use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufWriter, Read, Write},
    path::Path,
};

use anyhow::{Context, Result, anyhow, bail};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::io::Write as _,
    alignment::record::Flags,
    alignment::record::data::field::Tag,
    alignment::record_buf::data::field::Value,
    header::record::value::{Map, map::{Program, program::tag as pg_tag}},
};

/// Run the `flip` subcommand: read the RC-reference alignment, flip every record
/// back onto the original strand, and write the result as BAM.
pub fn run(input: Option<&Path>, output: Option<&Path>) -> Result<()> {
    let reader: Box<dyn Read> = match input {
        Some(p) if p != Path::new("-") => Box::new(
            File::open(p).with_context(|| format!("failed to open input {}", p.display()))?,
        ),
        _ => Box::new(io::stdin().lock()),
    };
    let mut reader = bam::io::Reader::new(reader);
    let mut header = reader.read_header().context("failed to read BAM header")?;

    // Contig lengths indexed by reference sequence id, and a name -> length map for SA.
    let mut lens: Vec<usize> = Vec::with_capacity(header.reference_sequences().len());
    let mut len_by_name: HashMap<String, usize> = HashMap::new();
    for (name, rs) in header.reference_sequences() {
        let l = rs.length().get();
        lens.push(l);
        len_by_name.insert(String::from_utf8_lossy(name).into_owned(), l);
    }

    append_pg(&mut header)?;

    let writer: Box<dyn Write> = match output {
        Some(p) if p != Path::new("-") => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("failed to create output {}", p.display()))?,
        )),
        _ => Box::new(BufWriter::new(io::stdout().lock())),
    };
    let mut writer = bam::io::Writer::new(writer);
    writer
        .write_header(&header)
        .context("failed to write BAM header")?;

    let mut warned = false;
    for (i, result) in reader.record_bufs(&header).enumerate() {
        let mut record = result.with_context(|| format!("failed to read record {i}"))?;
        transform_record(&mut record, &lens, &len_by_name, &mut warned)
            .with_context(|| format!("failed to transform record {i}"))?;
        writer
            .write_alignment_record(&header, &record)
            .with_context(|| format!("failed to write record {i}"))?;
    }

    writer.try_finish().context("failed to finalize BAM output")?;
    Ok(())
}

fn append_pg(header: &mut sam::Header) -> Result<()> {
    let cl = std::env::args().collect::<Vec<_>>().join(" ");
    let program = Map::<Program>::builder()
        .insert(pg_tag::NAME, "map-reval")
        .insert(pg_tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(pg_tag::COMMAND_LINE, cl)
        .build()
        .context("failed to build @PG record")?;
    header
        .programs_mut()
        .add("map-reval", program)
        .context("failed to append @PG record")?;
    Ok(())
}

/// Transform one record from the RC reference onto the original strand.
fn transform_record(
    rec: &mut sam::alignment::RecordBuf,
    lens: &[usize],
    len_by_name: &HashMap<String, usize>,
    warned: &mut bool,
) -> Result<()> {
    warn_unknown_tags(rec, warned);

    let flags = rec.flags();
    let paired = flags.is_segmented();
    let self_mapped = !flags.is_unmapped();
    let mate_mapped = paired && !flags.is_mate_unmapped();

    if self_mapped {
        // Coordinate flip using this record's own reference span.
        if let (Some(ref_id), Some(start)) = (rec.reference_sequence_id(), rec.alignment_start()) {
            let l = contig_len(lens, ref_id)?;
            let span = rec.cigar().alignment_span();
            *rec.alignment_start_mut() = Some(flip_pos(l, start.get(), span)?);
        }
        rec.cigar_mut().as_mut().reverse();
        rec.flags_mut().toggle(Flags::REVERSE_COMPLEMENTED);
        revcomp(rec.sequence_mut().as_mut());
        rec.quality_scores_mut().as_mut().reverse();

        // MD:Z -> reverse tokens, complement bases.
        if let Some(md) = get_string(rec, b"MD") {
            let new_md = reverse_md(&md);
            rec.data_mut()
                .insert(tag(b"MD"), Value::String(new_md.into()));
        }
        // SA:Z -> flip each supplementary alignment.
        if let Some(sa) = get_string(rec, b"SA") {
            let new_sa = transform_sa(&sa, len_by_name)?;
            rec.data_mut()
                .insert(tag(b"SA"), Value::String(new_sa.into_bytes().into()));
        }
    } else if mate_mapped {
        // Unmapped read placed at its mapped mate's coordinate: flip the placement
        // using the mate CIGAR so it stays co-located with the flipped mate.
        if let (Some(ref_id), Some(start), Some(mc)) =
            (rec.reference_sequence_id(), rec.alignment_start(), get_string(rec, b"MC"))
        {
            let l = contig_len(lens, ref_id)?;
            let span = cigar_ref_span(&mc)?;
            *rec.alignment_start_mut() = Some(flip_pos(l, start.get(), span)?);
        }
    }

    if paired {
        if mate_mapped {
            rec.flags_mut().toggle(Flags::MATE_REVERSE_COMPLEMENTED);

            let mate_ref = rec.mate_reference_sequence_id();
            let pnext = rec.mate_alignment_start();
            if let (Some(mate_ref), Some(pnext)) = (mate_ref, pnext) {
                if let Some(mc) = get_string(rec, b"MC") {
                    let l = contig_len(lens, mate_ref)?;
                    let span = cigar_ref_span(&mc)?;
                    *rec.mate_alignment_start_mut() = Some(flip_pos(l, pnext.get(), span)?);
                    let new_mc = reverse_cigar_str(&mc)?;
                    rec.data_mut()
                        .insert(tag(b"MC"), Value::String(new_mc.into_bytes().into()));
                } else if !*warned {
                    eprintln!(
                        "map-reval flip: warning: paired record lacks an MC tag; PNEXT left unchanged (further warnings suppressed)"
                    );
                    *warned = true;
                }
            }
        }
        *rec.template_length_mut() = -rec.template_length();
    }

    Ok(())
}

fn contig_len(lens: &[usize], id: usize) -> Result<usize> {
    lens.get(id)
        .copied()
        .ok_or_else(|| anyhow!("reference sequence id {id} out of range"))
}

fn get_string(rec: &sam::alignment::RecordBuf, key: &[u8; 2]) -> Option<Vec<u8>> {
    match rec.data().get(&tag(key)) {
        Some(Value::String(s)) => {
            let bytes: &[u8] = s.as_ref();
            Some(bytes.to_vec())
        }
        _ => None,
    }
}

fn tag(key: &[u8; 2]) -> Tag {
    Tag::from(*key)
}

// ---------------------------------------------------------------------------
// Pure helpers
// ---------------------------------------------------------------------------

/// New 1-based POS after a per-contig reverse complement: `L - pos - span + 2`.
fn flip_pos(l: usize, pos: usize, span: usize) -> Result<Position> {
    let end = pos
        .checked_add(span)
        .filter(|&e| e <= l + 1)
        .ok_or_else(|| anyhow!("alignment (pos={pos}, span={span}) runs past contig length {l}"))?;
    let new = l + 2 - end;
    Position::new(new).ok_or_else(|| anyhow!("computed a non-positive position"))
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', b'U' => b'A',
        b'R' => b'Y', b'Y' => b'R', b'S' => b'S', b'W' => b'W',
        b'K' => b'M', b'M' => b'K', b'B' => b'V', b'V' => b'B',
        b'D' => b'H', b'H' => b'D', b'N' => b'N',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c', b'u' => b'a',
        b'r' => b'y', b'y' => b'r', b's' => b's', b'w' => b'w',
        b'k' => b'm', b'm' => b'k', b'b' => b'v', b'v' => b'b',
        b'd' => b'h', b'h' => b'd', b'n' => b'n',
        other => other,
    }
}

fn revcomp(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement(*b);
    }
}

/// Sum of the reference-consuming CIGAR ops (M, D, N, =, X) in a CIGAR string.
fn cigar_ref_span(cigar: &[u8]) -> Result<usize> {
    let mut span = 0usize;
    for (len, op) in parse_cigar(cigar)? {
        if matches!(op, b'M' | b'D' | b'N' | b'=' | b'X') {
            span += len;
        }
    }
    Ok(span)
}

/// Reverse the operation order of a CIGAR string (e.g. `10S90M` -> `90M10S`).
fn reverse_cigar_str(cigar: &[u8]) -> Result<String> {
    if cigar == b"*" {
        return Ok("*".to_string());
    }
    let mut ops = parse_cigar(cigar)?;
    ops.reverse();
    let mut out = String::new();
    for (len, op) in ops {
        out.push_str(&len.to_string());
        out.push(op as char);
    }
    Ok(out)
}

fn parse_cigar(cigar: &[u8]) -> Result<Vec<(usize, u8)>> {
    let mut ops = Vec::new();
    let mut num: usize = 0;
    let mut seen_digit = false;
    for &c in cigar {
        if c.is_ascii_digit() {
            num = num * 10 + usize::from(c - b'0');
            seen_digit = true;
        } else if c.is_ascii_alphabetic() || c == b'=' {
            if !seen_digit {
                bail!("malformed CIGAR: operation without length");
            }
            ops.push((num, c));
            num = 0;
            seen_digit = false;
        } else {
            bail!("malformed CIGAR: unexpected byte {c:#x}");
        }
    }
    if seen_digit {
        bail!("malformed CIGAR: trailing length without operation");
    }
    Ok(ops)
}

/// Reverse an MD string: reverse token order, complement mismatch and deletion
/// bases, and reverse the base order within each deletion group.
/// e.g. `10A5^AC3` -> `3^GT5T10`.
fn reverse_md(md: &[u8]) -> Vec<u8> {
    enum Tok {
        Num(usize),
        Sub(u8),
        Del(Vec<u8>),
    }

    let mut toks = Vec::new();
    let mut i = 0;
    while i < md.len() {
        let c = md[i];
        if c.is_ascii_digit() {
            let mut n = 0usize;
            while i < md.len() && md[i].is_ascii_digit() {
                n = n * 10 + usize::from(md[i] - b'0');
                i += 1;
            }
            toks.push(Tok::Num(n));
        } else if c == b'^' {
            i += 1;
            let mut bases = Vec::new();
            while i < md.len() && md[i].is_ascii_alphabetic() {
                bases.push(md[i]);
                i += 1;
            }
            toks.push(Tok::Del(bases));
        } else if c.is_ascii_alphabetic() {
            toks.push(Tok::Sub(c));
            i += 1;
        } else {
            // Unexpected byte: pass through as a single-base substitution-like token.
            toks.push(Tok::Sub(c));
            i += 1;
        }
    }

    let mut out = Vec::new();
    for tok in toks.into_iter().rev() {
        match tok {
            Tok::Num(n) => out.extend_from_slice(n.to_string().as_bytes()),
            Tok::Sub(b) => out.push(complement(b)),
            Tok::Del(bases) => {
                out.push(b'^');
                for &b in bases.iter().rev() {
                    out.push(complement(b));
                }
            }
        }
    }
    out
}

/// Transform an SA:Z tag. Each `rname,pos,strand,CIGAR,mapQ,NM;` entry gets its
/// position flipped, CIGAR reversed, and strand toggled.
fn transform_sa(sa: &[u8], len_by_name: &HashMap<String, usize>) -> Result<String> {
    let sa = std::str::from_utf8(sa).context("SA tag is not valid UTF-8")?;
    let mut out = String::new();
    for entry in sa.split_terminator(';') {
        if entry.is_empty() {
            continue;
        }
        let mut f = entry.split(',');
        let rname = f.next().context("SA entry missing rname")?;
        let pos: usize = f
            .next()
            .context("SA entry missing pos")?
            .parse()
            .context("SA entry has invalid pos")?;
        let strand = f.next().context("SA entry missing strand")?;
        let cigar = f.next().context("SA entry missing CIGAR")?;
        let mapq = f.next().context("SA entry missing mapQ")?;
        let nm = f.next().context("SA entry missing NM")?;

        let l = *len_by_name
            .get(rname)
            .ok_or_else(|| anyhow!("SA references unknown contig {rname}"))?;
        let span = cigar_ref_span(cigar.as_bytes())?;
        let new_pos = flip_pos(l, pos, span)?;
        let new_strand = match strand {
            "+" => "-",
            "-" => "+",
            other => other,
        };
        let new_cigar = reverse_cigar_str(cigar.as_bytes())?;

        out.push_str(rname);
        out.push(',');
        out.push_str(&new_pos.get().to_string());
        out.push(',');
        out.push_str(new_strand);
        out.push(',');
        out.push_str(&new_cigar);
        out.push(',');
        out.push_str(mapq);
        out.push(',');
        out.push_str(nm);
        out.push(';');
    }
    Ok(out)
}

/// Tags that are either explicitly transformed or known to be orientation-invariant.
const KNOWN_TAGS: &[[u8; 2]] = &[
    *b"MD", *b"MC", *b"SA", // transformed
    *b"NM", *b"AS", *b"XS", *b"MQ", *b"ms", *b"md", *b"RG", *b"PG", *b"NH", *b"HI",
];

fn warn_unknown_tags(rec: &sam::alignment::RecordBuf, warned: &mut bool) {
    if *warned {
        return;
    }
    for key in rec.data().keys() {
        let bytes: [u8; 2] = key.into();
        if !KNOWN_TAGS.contains(&bytes) {
            eprintln!(
                "map-reval flip: warning: tag {}{} is not recognized; passing through unchanged (further warnings suppressed)",
                bytes[0] as char, bytes[1] as char
            );
            *warned = true;
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flip_pos_chr2() {
        // chr2 LN 242,193,529; forward POS 33,031,251 + 151bp span.
        let p = flip_pos(242_193_529, 33_031_251, 151).unwrap();
        assert_eq!(p.get(), 209_162_129);
    }

    #[test]
    fn flip_pos_is_involution() {
        let l = 1000;
        let p = flip_pos(l, 100, 50).unwrap();
        let back = flip_pos(l, p.get(), 50).unwrap();
        assert_eq!(back.get(), 100);
    }

    #[test]
    fn flip_pos_rejects_overrun() {
        assert!(flip_pos(100, 90, 20).is_err());
    }

    #[test]
    fn revcomp_basic() {
        let mut s = b"ACGTN".to_vec();
        revcomp(&mut s);
        assert_eq!(s, b"NACGT");
    }

    #[test]
    fn reverse_cigar_moves_clips() {
        assert_eq!(reverse_cigar_str(b"10S90M").unwrap(), "90M10S");
        assert_eq!(reverse_cigar_str(b"5H10M2D30M5S").unwrap(), "5S30M2D10M5H");
        assert_eq!(reverse_cigar_str(b"*").unwrap(), "*");
    }

    #[test]
    fn cigar_span_counts_reference_ops() {
        assert_eq!(cigar_ref_span(b"10M2D30M5S").unwrap(), 42);
        assert_eq!(cigar_ref_span(b"151M").unwrap(), 151);
    }

    #[test]
    fn reverse_md_example() {
        assert_eq!(reverse_md(b"10A5^AC3"), b"3^GT5T10");
    }

    #[test]
    fn reverse_md_is_involution() {
        let md = b"6G4C1^ATG12A0T5";
        let once = reverse_md(md);
        let twice = reverse_md(&once);
        assert_eq!(twice, md);
    }

    #[test]
    fn transform_sa_entry() {
        let mut lens = HashMap::new();
        lens.insert("chr1".to_string(), 1000usize);
        // pos 101, 10S90M (span 90) -> flip_pos(1000,101,90)=1000+2-191=811; cigar 90M10S; strand flip.
        let out = transform_sa(b"chr1,101,+,10S90M,60,2;", &lens).unwrap();
        assert_eq!(out, "chr1,811,-,90M10S,60,2;");
    }
}
