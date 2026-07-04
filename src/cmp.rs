use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};

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

    // Legacy tallies indexed by q = max(mapqA, mapqB).
    let mut n_reads = [0u64; 256];
    let mut n_wrong = [0u64; 256];
    let mut n_unmapped = [0u64; 256];

    // A-group indexed by mapQ_A (reads mapped in A); how they fare in B.
    let mut a_reads = [0u64; 256];
    let mut a_diff = [0u64; 256];
    let mut a_unmap = [0u64; 256];
    // B-group indexed by mapQ_B (reads mapped in B); how they fare in A.
    let mut b_reads = [0u64; 256];
    let mut b_diff = [0u64; 256];
    let mut b_unmap = [0u64; 256];

    let mut skipped_a = 0u64;
    let mut skipped_b = 0u64;
    let mut spliced = 0u64;
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
        spliced += u64::from(pa.spliced) + u64::from(pb.spliced);

        let qa = pa.mapq as usize;
        let qb = pb.mapq as usize;
        let q = qa.max(qb);

        let discordant = if pa.mapped && pb.mapped {
            let ratio = if pa.contig == pb.contig {
                reciprocal_overlap(pa.start, pa.end, pb.start, pb.end)
            } else {
                0.0
            };
            let wrong = ratio < min_overlap;
            n_reads[q] += 1;
            a_reads[qa] += 1;
            b_reads[qb] += 1;
            if wrong {
                n_wrong[q] += 1;
                a_diff[qa] += 1;
                b_diff[qb] += 1;
            }
            wrong
        } else if pa.mapped {
            // A mapped, B unmapped.
            n_reads[q] += 1;
            n_unmapped[q] += 1;
            a_reads[qa] += 1;
            a_unmap[qa] += 1;
            true
        } else if pb.mapped {
            // B mapped, A unmapped.
            n_reads[q] += 1;
            n_unmapped[q] += 1;
            b_reads[qb] += 1;
            b_unmap[qb] += 1;
            true
        } else {
            // Both unmapped: no placement to report; tallied on the U line.
            both_unmapped += 1;
            false
        };

        if emit_e && discordant {
            let name = String::from_utf8_lossy(rec_a.name().unwrap_or_default());
            writeln!(out, "E\t{name}\t{}\t{}", pa.fields(), pb.fields())
                .context("failed to write E line")?;
        }

        idx += 1;
    }

    for q in (0..256).rev() {
        if a_reads[q] > 0 || b_reads[q] > 0 || n_reads[q] > 0 {
            writeln!(
                out,
                "Q\t{q}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                a_reads[q], a_diff[q], a_unmap[q],
                b_reads[q], b_diff[q], b_unmap[q],
                n_reads[q], n_wrong[q], n_unmapped[q],
            )
            .context("failed to write Q line")?;
        }
    }
    writeln!(out, "U\t{both_unmapped}").context("failed to write U line")?;

    out.flush().context("failed to flush output")?;

    eprintln!(
        "map-reval cmp: compared {idx} primary pairs; skipped non-primary A={skipped_a} B={skipped_b}"
    );
    if spliced > 0 {
        eprintln!(
            "map-reval cmp: warning: {spliced} spliced primary alignment(s) with N ops; overlap is approximate (splice-aware comparison not yet implemented)"
        );
    }

    Ok(())
}

/// Reference-interval placement of a primary alignment on its contig.
struct Placement<'a> {
    mapped: bool,
    contig: Option<&'a str>,
    start: usize, // half-open interval [start, end) for overlap
    end: usize,
    rev: bool,
    mapq: u8,
    spliced: bool,
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
                spliced: false,
            };
        }
        let start = rec.alignment_start().map(|p| p.get()).unwrap_or(0);
        let span = rec.cigar().alignment_span();
        let contig = rec
            .reference_sequence_id()
            .and_then(|id| names.get(id))
            .map(String::as_str);
        let spliced = rec
            .cigar()
            .as_ref()
            .iter()
            .any(|op| op.kind() == sam::alignment::record::cigar::op::Kind::Skip);
        Placement {
            mapped: true,
            contig,
            start,
            end: start + span,
            rev: rec.flags().is_reverse_complemented(),
            mapq: rec.mapping_quality().map(|m| m.get()).unwrap_or(0),
            spliced,
        }
    }

    /// The five TAB-separated E-line fields for this side (`.` when unmapped).
    fn fields(&self) -> String {
        if !self.mapped {
            return ".\t.\t.\t.\t.".to_string();
        }
        let ctg = self.contig.unwrap_or(".");
        let strand = if self.rev { '-' } else { '+' };
        // Display 1-based inclusive end.
        format!("{ctg}\t{}\t{}\t{strand}\t{}", self.start, self.end - 1, self.mapq)
    }
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

/// Reciprocal overlap (intersection / union) of two half-open intervals on the
/// same contig. Callers pass 0-length/disjoint intervals for the no-overlap case.
fn reciprocal_overlap(a0: usize, a1: usize, b0: usize, b1: usize) -> f64 {
    let inter = a1.min(b1).saturating_sub(a0.max(b0));
    let union = (a1 - a0) + (b1 - b0) - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlap_identical() {
        assert_eq!(reciprocal_overlap(100, 200, 100, 200), 1.0);
    }

    #[test]
    fn overlap_half() {
        // [100,200) vs [150,250): inter=50, union=150 -> 1/3.
        let r = reciprocal_overlap(100, 200, 150, 250);
        assert!((r - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn overlap_disjoint() {
        assert_eq!(reciprocal_overlap(100, 200, 300, 400), 0.0);
    }

    #[test]
    fn overlap_containment() {
        // [100,200) contains [120,160): inter=40, union=100 -> 0.4.
        let r = reciprocal_overlap(100, 200, 120, 160);
        assert!((r - 0.4).abs() < 1e-9);
    }

    #[test]
    fn overlap_touching_is_zero() {
        assert_eq!(reciprocal_overlap(100, 200, 200, 300), 0.0);
    }
}
