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

    // Q: read-level concordance by reciprocal block-set overlap.
    //   *_a binned by mapQ_A, *_b by mapQ_B, *_m by max(mapQ_A, mapQ_B).
    let mut q_a = Group::new();
    let mut q_b = Group::new();
    let mut q_m = Group::new();
    // I: read-level intron-chain concordance (spliced reads only).
    let mut i_a = Group::new();
    let mut i_b = Group::new();
    let mut i_m = Group::new();
    // J: per-junction concordance.
    let mut j_a = Group::new();
    let mut j_b = Group::new();

    let mut skipped_a = 0u64;
    let mut skipped_b = 0u64;
    let mut spliced_a = 0u64;
    let mut spliced_b = 0u64;
    let mut both_unmapped = 0u64;

    let mut it_a = ra.record_bufs(&ha);
    let mut it_b = rb.record_bufs(&hb);
    let mut idx = 0u64;

    loop {
        let ca = next_primary(&mut it_a, &mut skipped_a).context("reading A")?;
        let cb = next_primary(&mut it_b, &mut skipped_b).context("reading B")?;

        let (rec_a, rec_b) = match (ca, cb) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => bail!(
                "primary count mismatch (A and B disagree after {idx} pairs); \
                 inputs must contain the same reads in the same order"
            ),
            (Some(a), Some(b)) => (a, b),
        };

        // Guard: the two streams must be the same read set in the same order.
        if rec_a.name() != rec_b.name()
            || rec_a.flags().is_first_segment() != rec_b.flags().is_first_segment()
            || rec_a.flags().is_last_segment() != rec_b.flags().is_last_segment()
        {
            bail!(
                "record mismatch at pair {idx}: A={:?} B={:?}; \
                 inputs must contain the same reads in the same order",
                String::from_utf8_lossy(rec_a.name().unwrap_or_default()),
                String::from_utf8_lossy(rec_b.name().unwrap_or_default()),
            );
        }

        let pa = Placement::of(&rec_a, &names_a);
        let pb = Placement::of(&rec_b, &names_b);

        let qa = pa.mapq as usize;
        let qb = pb.mapq as usize;
        let q = qa.max(qb);

        let same_contig = pa.mapped && pb.mapped && pa.contig == pb.contig;

        // ---- Q: reciprocal block-set overlap ----
        let discordant = if pa.mapped && pb.mapped {
            let ratio = if same_contig {
                reciprocal_overlap(&pa.blocks, &pb.blocks)
            } else {
                0.0
            };
            let wrong = ratio < min_overlap;
            q_m.add(q, wrong, false);
            q_a.add(qa, wrong, false);
            q_b.add(qb, wrong, false);
            wrong
        } else if pa.mapped {
            q_m.add(q, false, true);
            q_a.add(qa, false, true);
            true
        } else if pb.mapped {
            q_m.add(q, false, true);
            q_b.add(qb, false, true);
            true
        } else {
            both_unmapped += 1;
            false
        };

        // ---- I: intron-chain equality (spliced reads only) ----
        let a_spliced = !pa.junctions.is_empty();
        let b_spliced = !pb.junctions.is_empty();
        spliced_a += u64::from(a_spliced);
        spliced_b += u64::from(b_spliced);
        if a_spliced || b_spliced {
            let chain_ok = same_contig && pa.junctions == pb.junctions;
            if a_spliced {
                i_a.add(qa, pb.mapped && !chain_ok, !pb.mapped);
            }
            if b_spliced {
                i_b.add(qb, pa.mapped && !chain_ok, !pa.mapped);
            }
            let one_unmapped = !pa.mapped || !pb.mapped;
            i_m.add(q, !one_unmapped && !chain_ok, one_unmapped);
        }

        // ---- J: per-junction exact match (linear via sorted-list merge) ----
        let shared = if same_contig {
            shared_count(&pa.junctions, &pb.junctions) as u64
        } else {
            0
        };
        let na = pa.junctions.len() as u64;
        if na > 0 {
            // na > 0 implies A is mapped; a junction is `diff` when B is mapped
            // but lacks an exact match, `unmap` when B is unmapped.
            if pb.mapped {
                j_a.add_n(qa, na, na - shared, 0);
            } else {
                j_a.add_n(qa, na, 0, na);
            }
        }
        let nb = pb.junctions.len() as u64;
        if nb > 0 {
            if pa.mapped {
                j_b.add_n(qb, nb, nb - shared, 0);
            } else {
                j_b.add_n(qb, nb, 0, nb);
            }
        }

        if emit_e && discordant {
            let name = String::from_utf8_lossy(rec_a.name().unwrap_or_default());
            writeln!(out, "E\t{name}\t{}\t{}", pa.fields(), pb.fields())
                .context("failed to write E line")?;
        }

        idx += 1;
    }

    write_grouped(&mut out, "Q", &q_a, &q_b, Some(&q_m))?;
    write_grouped(&mut out, "I", &i_a, &i_b, Some(&i_m))?;
    write_grouped(&mut out, "J", &j_a, &j_b, None)?;
    writeln!(out, "U\t{both_unmapped}").context("failed to write U line")?;

    out.flush().context("failed to flush output")?;

    eprintln!(
        "map-reval cmp: compared {idx} primary pairs (spliced A={spliced_a} B={spliced_b}); \
         skipped non-primary A={skipped_a} B={skipped_b}"
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

    fn add_n(&mut self, q: usize, reads: u64, diff: u64, unmap: u64) {
        self.reads[q] += reads;
        self.diff[q] += diff;
        self.unmap[q] += unmap;
    }
}

/// Print grouped per-mapQ rows (high→low) for a line type. `m` is the
/// max(mapQ)-binned trio, appended when present (Q, I) and omitted for J.
fn write_grouped(
    out: &mut dyn Write,
    tag: &str,
    a: &Group,
    b: &Group,
    m: Option<&Group>,
) -> Result<()> {
    for q in (0..256).rev() {
        let show = a.reads[q] > 0 || b.reads[q] > 0 || m.is_some_and(|m| m.reads[q] > 0);
        if !show {
            continue;
        }
        write!(
            out,
            "{tag}\t{q}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.reads[q], a.diff[q], a.unmap[q], b.reads[q], b.diff[q], b.unmap[q],
        )
        .context("failed to write summary line")?;
        if let Some(m) = m {
            write!(out, "\t{}\t{}\t{}", m.reads[q], m.diff[q], m.unmap[q])
                .context("failed to write summary line")?;
        }
        writeln!(out).context("failed to write summary line")?;
    }
    Ok(())
}

/// A list of half-open reference intervals (exon blocks or intron junctions).
type Intervals = Vec<(usize, usize)>;

/// Reference placement of a primary alignment: exon blocks and intron junctions.
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

    /// The five TAB-separated E-line fields for this side (`.` when unmapped).
    fn fields(&self) -> String {
        if !self.mapped {
            return ".\t.\t.\t.\t.".to_string();
        }
        let ctg = self.contig.unwrap_or(".");
        let strand = if self.rev { '-' } else { '+' };
        // Display 1-based inclusive end of the outer extent.
        format!("{ctg}\t{}\t{}\t{strand}\t{}", self.start, self.end - 1, self.mapq)
    }
}

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

fn next_primary<I>(it: &mut I, skipped: &mut u64) -> Result<Option<RecordBuf>>
where
    I: Iterator<Item = io::Result<RecordBuf>>,
{
    for r in it.by_ref() {
        let rec = r?;
        let f = rec.flags();
        if f.is_secondary() || f.is_supplementary() {
            *skipped += 1;
            continue;
        }
        return Ok(Some(rec));
    }
    Ok(None)
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
    fn overlap_two_block_partial() {
        // A exons [100,150)+[350,400) (len 100); B exons [100,150)+[360,400) (len 90).
        // intersection = 50 + 40 = 90; union = 100 + 90 - 90 = 100 -> 0.9.
        let a = [(100, 150), (350, 400)];
        let b = [(100, 150), (360, 400)];
        let r = reciprocal_overlap(&a, &b);
        assert!((r - 0.9).abs() < 1e-9);
    }
}
