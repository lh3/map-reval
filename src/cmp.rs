use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, record::cigar::op::Kind},
};

/// Compare two BAMs derived from the same reads in the same order, on the same
/// coordinate system (typically a forward-reference alignment vs a `flip`ped
/// reverse-complement alignment). Emits a TAB-delimited, line-typed report.
pub fn run(
    a: &Path,
    b: &Path,
    min_overlap: f64,
    output: Option<&Path>,
    emit_e: bool,
    min_mapq: u8,
) -> Result<()> {
    let mut ra = bam::io::Reader::new(
        File::open(a).with_context(|| format!("failed to open {}", a.display()))?,
    );
    let ha = ra.read_header().context("failed to read header of A")?;
    let names_a = contig_names(&ha);

    let mut rb = bam::io::Reader::new(
        File::open(b).with_context(|| format!("failed to open {}", b.display()))?,
    );
    let hb = rb.read_header().context("failed to read header of B")?;
    let names_b = contig_names(&hb);

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("failed to create {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    let mut cx = Counters::new();
    let mut secondary_a = 0u64;
    let mut secondary_b = 0u64;
    let mut supp_a = 0u64;
    let mut supp_b = 0u64;

    let mut it_a = ra.record_bufs(&ha);
    let mut it_b = rb.record_bufs(&hb);
    let mut pending_a: Option<RecordBuf> = None;
    let mut pending_b: Option<RecordBuf> = None;
    let mut idx = 0u64;

    loop {
        let ga = next_group(&mut it_a, &mut pending_a).context("reading A")?;
        let gb = next_group(&mut it_b, &mut pending_b).context("reading B")?;

        let (group_a, group_b) = match (ga, gb) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => bail!(
                "read-group count mismatch (A and B disagree after {idx} pairs); \
                 inputs must contain the same reads in the same order"
            ),
            (Some(a), Some(b)) => (a, b),
        };

        if group_name(&group_a) != group_name(&group_b) {
            bail!(
                "record mismatch at pair {idx}: A={:?} B={:?}; \
                 inputs must contain the same reads in the same order",
                String::from_utf8_lossy(group_name(&group_a).unwrap_or_default()),
                String::from_utf8_lossy(group_name(&group_b).unwrap_or_default()),
            );
        }

        let segs_a = split_group(&group_a, &names_a, &mut secondary_a, &mut supp_a);
        let segs_b = split_group(&group_b, &names_b, &mut secondary_b, &mut supp_b);

        for key in 0..segs_a.len() {
            match (&segs_a[key], &segs_b[key]) {
                (None, None) => continue,
                (Some(_), None) | (None, Some(_)) => bail!(
                    "segment mismatch at pair {idx}: read present in only one file; \
                     inputs must contain the same reads in the same order"
                ),
                (Some(sa), Some(sb)) => {
                    let discordant = cx.compare(sa, sb, min_overlap);

                    let has_junctions =
                        !sa.primary.junctions.is_empty() || !sb.primary.junctions.is_empty();
                    if emit_e
                        && sa.primary.mapq.max(sb.primary.mapq) >= min_mapq
                        && (discordant || has_junctions)
                    {
                        let name = format!("{}{}", group_name_str(&group_a), seg_suffix(key));
                        if discordant {
                            writeln!(
                                out,
                                "E\t{name}\t{}\t{}",
                                sa.primary.fields(),
                                sb.primary.fields()
                            )
                            .context("failed to write E line")?;
                        }
                        let same_contig = sa.primary.mapped
                            && sb.primary.mapped
                            && sa.primary.contig == sb.primary.contig;
                        write_f_lines(&mut out, &name, &sa.primary, &sb.primary, same_contig)?;
                    }
                    idx += 1;
                }
            }
        }
    }

    write_grouped(&mut out, "Q", &cx.q_a, &cx.q_b)?;
    write_grouped(&mut out, "I", &cx.i_a, &cx.i_b)?;
    write_j(&mut out, &cx.j_a, &cx.j_b)?;
    writeln!(out, "U\t{}", cx.both_unmapped).context("failed to write U line")?;

    out.flush().context("failed to flush output")?;

    eprintln!(
        "map-reval cmp: compared {idx} primary pairs (spliced A={} B={}); \
         supplementary A={supp_a} B={supp_b}; secondary A={secondary_a} B={secondary_b}",
        cx.spliced_a, cx.spliced_b
    );

    Ok(())
}

/// Per-mapQ tallies for one column group: total, discordant, and unmapped-in-other.
struct Group {
    reads: [u64; 256],
    diff: [u64; 256],
    unmap: [u64; 256],
}

impl Group {
    fn new() -> Self {
        Self {
            reads: [0; 256],
            diff: [0; 256],
            unmap: [0; 256],
        }
    }

    fn add(&mut self, q: usize, diff: bool, unmap: bool) {
        self.reads[q] += 1;
        if diff {
            self.diff[q] += 1;
        }
        if unmap {
            self.unmap[q] += 1;
        }
    }
}

/// Per-mapQ junction tallies. `gone` = other read mapped but no overlapping
/// junction; `shifted` = overlaps a junction but not an exact match; `unmap` =
/// other read unmapped. Exact matches = `at - gone - shifted - unmap`.
struct JGroup {
    at: [u64; 256],
    gone: [u64; 256],
    shifted: [u64; 256],
    unmap: [u64; 256],
}

impl JGroup {
    fn new() -> Self {
        Self {
            at: [0; 256],
            gone: [0; 256],
            shifted: [0; 256],
            unmap: [0; 256],
        }
    }

    fn add_n(&mut self, q: usize, at: u64, gone: u64, shifted: u64, unmap: u64) {
        self.at[q] += at;
        self.gone[q] += gone;
        self.shifted[q] += shifted;
        self.unmap[q] += unmap;
    }
}

/// Print grouped per-mapQ rows (high→low) for a line type: an A-side and a
/// B-side triple, binned by each file's own mapQ.
fn write_grouped(out: &mut dyn Write, tag: &str, a: &Group, b: &Group) -> Result<()> {
    for q in (0..256).rev() {
        if a.reads[q] == 0 && b.reads[q] == 0 {
            continue;
        }
        writeln!(
            out,
            "{tag}\t{q}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.reads[q], a.diff[q], a.unmap[q], b.reads[q], b.diff[q], b.unmap[q],
        )
        .context("failed to write summary line")?;
    }
    Ok(())
}

/// Emit an F line per non-exact junction on both sides (only under `-e`).
/// The focus junction fills its own side; the other side shows the largest-
/// overlap partner (shifted) or `.` (gone/unmapped).
fn write_f_lines(
    out: &mut dyn Write,
    name: &str,
    pa: &Placement,
    pb: &Placement,
    same_contig: bool,
) -> Result<()> {
    let mut lines: Vec<String> = Vec::new();
    for &j in &pa.junctions {
        let class = if pb.mapped && same_contig {
            classify(j, &pb.junctions)
        } else {
            JClass::Gone
        };
        match class {
            JClass::Exact => {}
            JClass::Shifted(o) => {
                lines.push(format!("F\t{name}\t{}\t{}", pa.junc_fields(j), pb.junc_fields(o)))
            }
            JClass::Gone => lines.push(format!("F\t{name}\t{}\t{DOT5}", pa.junc_fields(j))),
        }
    }
    for &j in &pb.junctions {
        let class = if pa.mapped && same_contig {
            classify(j, &pa.junctions)
        } else {
            JClass::Gone
        };
        match class {
            JClass::Exact => {}
            JClass::Shifted(o) => {
                lines.push(format!("F\t{name}\t{}\t{}", pa.junc_fields(o), pb.junc_fields(j)))
            }
            JClass::Gone => lines.push(format!("F\t{name}\t{DOT5}\t{}", pb.junc_fields(j))),
        }
    }
    // A symmetric shifted pair yields two identical lines; emit each unique line once.
    let mut seen = std::collections::HashSet::new();
    for line in &lines {
        if seen.insert(line.as_str()) {
            writeln!(out, "{line}").context("failed to write F line")?;
        }
    }
    Ok(())
}

/// Print the J (per-junction) lines: 8 data columns, no max-binned trio.
fn write_j(out: &mut dyn Write, a: &JGroup, b: &JGroup) -> Result<()> {
    for q in (0..256).rev() {
        if a.at[q] == 0 && b.at[q] == 0 {
            continue;
        }
        writeln!(
            out,
            "J\t{q}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.at[q], a.shifted[q], a.gone[q], a.unmap[q],
            b.at[q], b.shifted[q], b.gone[q], b.unmap[q],
        )
        .context("failed to write J line")?;
    }
    Ok(())
}

/// A list of half-open reference intervals (exon blocks or intron junctions).
type Intervals = Vec<(usize, usize)>;

/// Reference placement of one alignment: exon blocks and intron junctions.
#[derive(Clone)]
struct Placement<'a> {
    mapped: bool,
    contig: Option<&'a str>,
    start: usize, // outer extent, for the E-line display only
    end: usize,
    rev: bool,
    mapq: u8,
    blocks: Intervals,    // exon intervals [start, end), split at N
    junctions: Intervals, // intron intervals between blocks
}

impl<'a> Placement<'a> {
    fn of(rec: &RecordBuf, names: &'a [String]) -> Self {
        if rec.flags().is_unmapped() {
            return Placement {
                mapped: false,
                contig: None,
                start: 0,
                end: 0,
                rev: false,
                mapq: 0,
                blocks: Vec::new(),
                junctions: Vec::new(),
            };
        }
        let start = rec.alignment_start().map(|p| p.get()).unwrap_or(1);
        let (blocks, junctions) =
            ref_blocks(start, rec.cigar().as_ref().iter().map(|op| (op.kind(), op.len())));
        let contig = rec
            .reference_sequence_id()
            .and_then(|id| names.get(id))
            .map(String::as_str);
        let end = blocks.last().map(|&(_, e)| e).unwrap_or(start);
        Placement {
            mapped: true,
            contig,
            start,
            end,
            rev: rec.flags().is_reverse_complemented(),
            // Raw MAPQ byte. noodles decodes the SAM "missing" sentinel (255) to
            // None; aligners like STAR use 255 for unique mappings, so recover it
            // as 255 rather than collapsing to 0 (unmapped reads use 0 above).
            mapq: rec.mapping_quality().map(|m| m.get()).unwrap_or(255),
            blocks,
            junctions,
        }
    }

    /// The five TAB-separated E-line fields for this side (`.` when unmapped):
    /// the outer extent as a 0-based half-open (BED) interval.
    fn fields(&self) -> String {
        if !self.mapped {
            return DOT5.to_string();
        }
        bed_fields(self.contig.unwrap_or("."), self.start, self.end, self.rev, self.mapq)
    }

    /// F-line fields for one junction of this read (BED coordinates).
    fn junc_fields(&self, j: (usize, usize)) -> String {
        bed_fields(self.contig.unwrap_or("."), j.0, j.1, self.rev, self.mapq)
    }
}

/// Five TAB-separated fields for one alignment side: a 0-based half-open (BED)
/// interval plus its read's contig/strand/mapQ. `lo`/`hi` are 1-based `[lo, hi)`.
fn bed_fields(ctg: &str, lo: usize, hi: usize, rev: bool, mapq: u8) -> String {
    let strand = if rev { '-' } else { '+' };
    format!("{ctg}\t{}\t{}\t{strand}\t{}", lo - 1, hi - 1, mapq)
}

/// The unmapped placeholder for a five-field alignment side.
const DOT5: &str = ".\t.\t.\t.\t.";

/// Split a CIGAR into exon blocks and intron junctions given a 1-based start.
/// `N` (Skip) ends the current exon and opens an intron; other reference-
/// consuming ops (M/D/=/X) extend the current exon; read-only ops are ignored.
fn ref_blocks<I>(start: usize, ops: I) -> (Intervals, Intervals)
where
    I: IntoIterator<Item = (Kind, usize)>,
{
    let mut blocks = Vec::new();
    let mut junctions = Vec::new();
    let mut pos = start;
    let mut block_start = start;
    for (kind, len) in ops {
        if kind == Kind::Skip {
            blocks.push((block_start, pos));
            junctions.push((pos, pos + len));
            pos += len;
            block_start = pos;
        } else if kind.consumes_reference() {
            pos += len;
        }
    }
    blocks.push((block_start, pos));
    (blocks, junctions)
}

fn contig_names(header: &sam::Header) -> Vec<String> {
    header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k).into_owned())
        .collect()
}

/// All per-mapQ tallies accumulated over the run.
struct Counters {
    q_a: Group,
    q_b: Group,
    i_a: Group,
    i_b: Group,
    j_a: JGroup,
    j_b: JGroup,
    spliced_a: u64,
    spliced_b: u64,
    both_unmapped: u64,
}

impl Counters {
    fn new() -> Self {
        Self {
            q_a: Group::new(),
            q_b: Group::new(),
            i_a: Group::new(),
            i_b: Group::new(),
            j_a: JGroup::new(),
            j_b: JGroup::new(),
            spliced_a: 0,
            spliced_b: 0,
            both_unmapped: 0,
        }
    }

    /// Compare one read (a single segment) and update all tallies; returns
    /// whether the pair is placement-discordant (drives E emission). `sa`/`sb`
    /// carry the primary placement plus every mapped alignment (primary +
    /// supplementary).
    fn compare(&mut self, sa: &Seg, sb: &Seg, min_overlap: f64) -> bool {
        let pa = &sa.primary;
        let pb = &sb.primary;
        let qa = pa.mapq as usize;
        let qb = pb.mapq as usize;
        let same_contig = pa.mapped && pb.mapped && pa.contig == pb.contig;

        // ---- Q: primary concordant if it overlaps ANY alignment on the other side ----
        let discordant = if pa.mapped && pb.mapped {
            let a_conc = sb.alns.iter().any(|o| {
                pa.contig == o.contig && reciprocal_overlap(&pa.blocks, &o.blocks) >= min_overlap
            });
            let b_conc = sa.alns.iter().any(|o| {
                pb.contig == o.contig && reciprocal_overlap(&pb.blocks, &o.blocks) >= min_overlap
            });
            self.q_a.add(qa, !a_conc, false);
            self.q_b.add(qb, !b_conc, false);
            !a_conc || !b_conc
        } else if pa.mapped {
            self.q_a.add(qa, false, true);
            true
        } else if pb.mapped {
            self.q_b.add(qb, false, true);
            true
        } else {
            self.both_unmapped += 1;
            false
        };

        // ---- I: primary chain concordant if it matches ANY alignment's chain ----
        let a_spliced = !pa.junctions.is_empty();
        let b_spliced = !pb.junctions.is_empty();
        self.spliced_a += u64::from(a_spliced);
        self.spliced_b += u64::from(b_spliced);
        if a_spliced || b_spliced {
            let a_chain_conc = sb
                .alns
                .iter()
                .any(|o| o.contig == pa.contig && o.junctions == pa.junctions);
            let b_chain_conc = sa
                .alns
                .iter()
                .any(|o| o.contig == pb.contig && o.junctions == pb.junctions);
            if a_spliced {
                self.i_a.add(qa, pb.mapped && !a_chain_conc, !pb.mapped);
            }
            if b_spliced {
                self.i_b.add(qb, pa.mapped && !b_chain_conc, !pa.mapped);
            }
        }

        // ---- J: per-junction match (primary-vs-primary, unchanged) ----
        let shared = if same_contig {
            shared_count(&pa.junctions, &pb.junctions) as u64
        } else {
            0
        };
        let na = pa.junctions.len() as u64;
        if na > 0 {
            if !pb.mapped {
                self.j_a.add_n(qa, na, 0, 0, na);
            } else if !same_contig {
                self.j_a.add_n(qa, na, na, 0, 0);
            } else {
                let ov = overlap_count(&pa.junctions, &pb.junctions) as u64;
                self.j_a.add_n(qa, na, na - ov, ov - shared, 0);
            }
        }
        let nb = pb.junctions.len() as u64;
        if nb > 0 {
            if !pa.mapped {
                self.j_b.add_n(qb, nb, 0, 0, nb);
            } else if !same_contig {
                self.j_b.add_n(qb, nb, nb, 0, 0);
            } else {
                let ov = overlap_count(&pb.junctions, &pa.junctions) as u64;
                self.j_b.add_n(qb, nb, nb - ov, ov - shared, 0);
            }
        }

        discordant
    }
}

/// One read segment: its primary placement plus every mapped alignment
/// (primary if mapped + supplementaries) for supplementary-aware concordance.
struct Seg<'a> {
    primary: Placement<'a>,
    alns: Vec<Placement<'a>>,
}

/// Read all consecutive records sharing a query name (contiguous under
/// `GO:query`), carrying the first record of the next group in `pending`.
fn next_group<I>(it: &mut I, pending: &mut Option<RecordBuf>) -> Result<Option<Vec<RecordBuf>>>
where
    I: Iterator<Item = io::Result<RecordBuf>>,
{
    let first = match pending.take() {
        Some(r) => r,
        None => match it.next() {
            Some(r) => r?,
            None => return Ok(None),
        },
    };
    let name: Option<Vec<u8>> = rec_name(&first).map(<[u8]>::to_vec);
    let mut group = vec![first];
    for r in it.by_ref() {
        let rec = r?;
        if rec_name(&rec) == name.as_deref() {
            group.push(rec);
        } else {
            *pending = Some(rec);
            break;
        }
    }
    Ok(Some(group))
}

fn rec_name(rec: &RecordBuf) -> Option<&[u8]> {
    rec.name().map(AsRef::as_ref)
}

fn group_name(group: &[RecordBuf]) -> Option<&[u8]> {
    group.first().and_then(rec_name)
}

fn group_name_str(group: &[RecordBuf]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(group_name(group).unwrap_or_default())
}

/// Segment key: 1 = first (0x40), 2 = last (0x80), 0 = unpaired.
fn seg_key(flags: sam::alignment::record::Flags) -> usize {
    if flags.is_first_segment() {
        1
    } else if flags.is_last_segment() {
        2
    } else {
        0
    }
}

fn seg_suffix(key: usize) -> &'static str {
    match key {
        1 => "/1",
        2 => "/2",
        _ => "",
    }
}

/// Split a name group into per-segment primary + mapped-alignment sets, skipping
/// secondary alignments (counted). `supp` accumulates supplementary counts.
fn split_group<'a>(
    group: &[RecordBuf],
    names: &'a [String],
    secondary: &mut u64,
    supp: &mut u64,
) -> [Option<Seg<'a>>; 3] {
    let mut segs: [Option<Seg>; 3] = [None, None, None];
    // Pass 1: primaries (records may be in any order within the name group).
    for rec in group {
        let f = rec.flags();
        if f.is_secondary() {
            *secondary += 1;
        } else if !f.is_supplementary() {
            let p = Placement::of(rec, names);
            let alns = if p.mapped { vec![p.clone()] } else { Vec::new() };
            segs[seg_key(f)] = Some(Seg { primary: p, alns });
        }
    }
    // Pass 2: attach mapped supplementaries to their segment's alignment set.
    for rec in group {
        let f = rec.flags();
        if f.is_secondary() || !f.is_supplementary() {
            continue;
        }
        *supp += 1;
        let p = Placement::of(rec, names);
        if !p.mapped {
            continue;
        }
        if let Some(seg) = segs[seg_key(f)].as_mut() {
            seg.alns.push(p);
        }
    }
    segs
}

/// Reciprocal overlap (intersection / union) of two exon-block sets on the same
/// contig. Each set is sorted and internally disjoint (natural CIGAR order).
fn reciprocal_overlap(a: &[(usize, usize)], b: &[(usize, usize)]) -> f64 {
    let sum_a: usize = a.iter().map(|&(s, e)| e - s).sum();
    let sum_b: usize = b.iter().map(|&(s, e)| e - s).sum();
    let mut inter = 0usize;
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let lo = a[i].0.max(b[j].0);
        let hi = a[i].1.min(b[j].1);
        if hi > lo {
            inter += hi - lo;
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    let union = sum_a + sum_b - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

/// Number of intervals present in both lists. Both must be sorted, strictly
/// ascending, and internally disjoint (as `ref_blocks` produces junctions).
fn shared_count(a: &[(usize, usize)], b: &[(usize, usize)]) -> usize {
    use std::cmp::Ordering::{Equal, Greater, Less};
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Less => i += 1,
            Greater => j += 1,
            Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// Classification of one junction against another alignment's junction set.
enum JClass {
    Exact,
    Shifted((usize, usize)),
    Gone,
}

/// Classify junction `j` against `others` (sorted, strictly ascending, disjoint):
/// exact-coordinate match, else the largest-overlap partner (shifted), else gone.
fn classify(j: (usize, usize), others: &[(usize, usize)]) -> JClass {
    let mut k = others.partition_point(|&(_, e)| e <= j.0); // first that can overlap
    let (mut best, mut best_ov) = (None, 0usize);
    while k < others.len() && others[k].0 < j.1 {
        let o = others[k];
        if o == j {
            return JClass::Exact;
        }
        let ov = j.1.min(o.1) - j.0.max(o.0);
        if ov > best_ov {
            best_ov = ov;
            best = Some(o);
        }
        k += 1;
    }
    best.map_or(JClass::Gone, JClass::Shifted)
}

/// Count intervals in `a` that intersect at least one interval in `b` (both
/// sorted, strictly ascending, disjoint). Each `a` interval is counted once.
fn overlap_count(a: &[(usize, usize)], b: &[(usize, usize)]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        if a[i].1 <= b[j].0 {
            i += 1; // a[i] entirely before b[j]
        } else if b[j].1 <= a[i].0 {
            j += 1; // b[j] entirely before a[i]
        } else {
            n += 1; // overlap: count a[i] once, advance a
            i += 1;
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ops(spec: &[(Kind, usize)]) -> Vec<(Kind, usize)> {
        spec.to_vec()
    }

    #[test]
    fn blocks_unspliced() {
        let (b, j) = ref_blocks(100, ops(&[(Kind::Match, 100)]));
        assert_eq!(b, vec![(100, 200)]);
        assert!(j.is_empty());
    }

    #[test]
    fn blocks_single_junction() {
        // 50M200N50M from 100 -> exons [100,150),[350,400); intron [150,350).
        let (b, j) = ref_blocks(
            100,
            ops(&[(Kind::Match, 50), (Kind::Skip, 200), (Kind::Match, 50)]),
        );
        assert_eq!(b, vec![(100, 150), (350, 400)]);
        assert_eq!(j, vec![(150, 350)]);
    }

    #[test]
    fn blocks_deletion_stays_in_exon() {
        // 10M5D10M2N10M from 100: exon [100,125), intron [125,127), exon [127,137).
        let (b, j) = ref_blocks(
            100,
            ops(&[
                (Kind::Match, 10),
                (Kind::Deletion, 5),
                (Kind::Match, 10),
                (Kind::Skip, 2),
                (Kind::Match, 10),
            ]),
        );
        assert_eq!(b, vec![(100, 125), (127, 137)]);
        assert_eq!(j, vec![(125, 127)]);
    }

    #[test]
    fn blocks_leading_softclip_ignored() {
        let (b, j) = ref_blocks(100, ops(&[(Kind::SoftClip, 5), (Kind::Match, 100)]));
        assert_eq!(b, vec![(100, 200)]);
        assert!(j.is_empty());
    }

    #[test]
    fn blocks_multi_junction() {
        // 31M118N45M143N45M143N29M from 1000 (three junctions).
        let (b, j) = ref_blocks(
            1000,
            ops(&[
                (Kind::Match, 31),
                (Kind::Skip, 118),
                (Kind::Match, 45),
                (Kind::Skip, 143),
                (Kind::Match, 45),
                (Kind::Skip, 143),
                (Kind::Match, 29),
            ]),
        );
        assert_eq!(j.len(), 3);
        assert_eq!(j[0], (1031, 1149));
        assert_eq!(b.len(), 4);
        assert_eq!(b[0], (1000, 1031));
        assert_eq!(b.last().unwrap().1, 1000 + 31 + 118 + 45 + 143 + 45 + 143 + 29);
    }

    // reciprocal_overlap over block sets; single-interval cases match the old scalar tests.
    #[test]
    fn overlap_identical() {
        assert_eq!(reciprocal_overlap(&[(100, 200)], &[(100, 200)]), 1.0);
    }

    #[test]
    fn overlap_half() {
        let r = reciprocal_overlap(&[(100, 200)], &[(150, 250)]);
        assert!((r - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn overlap_disjoint() {
        assert_eq!(reciprocal_overlap(&[(100, 200)], &[(300, 400)]), 0.0);
    }

    #[test]
    fn overlap_containment() {
        let r = reciprocal_overlap(&[(100, 200)], &[(120, 160)]);
        assert!((r - 0.4).abs() < 1e-9);
    }

    #[test]
    fn overlap_touching_is_zero() {
        assert_eq!(reciprocal_overlap(&[(100, 200)], &[(200, 300)]), 0.0);
    }

    #[test]
    fn overlap_two_block_identical() {
        let s = [(100, 150), (350, 400)];
        assert_eq!(reciprocal_overlap(&s, &s), 1.0);
    }

    #[test]
    fn shared_count_cases() {
        // disjoint
        assert_eq!(shared_count(&[(10, 20)], &[(30, 40)]), 0);
        // fully shared
        let s = [(10, 20), (30, 40), (50, 60)];
        assert_eq!(shared_count(&s, &s), 3);
        // partial with interleaving: shared (10,20) and (50,60); (30,40)/(35,45) differ.
        let a = [(10, 20), (30, 40), (50, 60)];
        let b = [(10, 20), (35, 45), (50, 60)];
        assert_eq!(shared_count(&a, &b), 2);
        // same start, different end must not match.
        assert_eq!(shared_count(&[(10, 20)], &[(10, 25)]), 0);
        // empty lists
        assert_eq!(shared_count(&[], &[(10, 20)]), 0);
    }

    #[test]
    fn classify_cases() {
        use JClass::*;
        assert!(matches!(classify((100, 200), &[(100, 200)]), Exact));
        assert!(matches!(classify((100, 200), &[(150, 250)]), Shifted((150, 250))));
        // multi-overlap picks the largest overlap: (150,250)∩=50 > (90,120)∩=20.
        assert!(matches!(
            classify((100, 200), &[(90, 120), (150, 250)]),
            Shifted((150, 250))
        ));
        assert!(matches!(classify((100, 200), &[(300, 400)]), Gone));
        assert!(matches!(classify((100, 200), &[]), Gone));
        // exact wins even when another junction also overlaps.
        assert!(matches!(classify((100, 200), &[(100, 200), (150, 250)]), Exact));
    }

    #[test]
    fn overlap_count_cases() {
        // disjoint
        assert_eq!(overlap_count(&[(10, 20)], &[(30, 40)]), 0);
        // touching is half-open -> no overlap
        assert_eq!(overlap_count(&[(10, 20)], &[(20, 30)]), 0);
        // one `a` spanning two `b` counts the single `a` once
        assert_eq!(overlap_count(&[(10, 100)], &[(20, 30), (40, 50)]), 1);
        // two disjoint `a` both overlapping one `b` counts both
        assert_eq!(overlap_count(&[(10, 22), (24, 40)], &[(20, 30)]), 2);
        // containment
        assert_eq!(overlap_count(&[(20, 30)], &[(10, 100)]), 1);
        // exact match is also an overlap
        assert_eq!(overlap_count(&[(10, 20)], &[(10, 20)]), 1);
        // empty
        assert_eq!(overlap_count(&[], &[(10, 20)]), 0);
    }

    #[test]
    fn overlap_two_block_partial() {
        // A exons [100,150)+[350,400) (len 100); B exons [100,150)+[360,400) (len 90).
        // intersection = 50 + 40 = 90; union = 100 + 90 - 90 = 100 -> 0.9.
        let a = [(100, 150), (350, 400)];
        let b = [(100, 150), (360, 400)];
        let r = reciprocal_overlap(&a, &b);
        assert!((r - 0.9).abs() < 1e-9);
    }
}
