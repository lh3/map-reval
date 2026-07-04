mod cmp;
mod flip;

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

/// Compare read alignments against a reference and its reverse-complement.
#[derive(Parser)]
#[command(name = "map-reval", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Flip an alignment made against the reverse-complement reference back onto
    /// the original strand, emitting a BAM comparable to a forward-reference alignment.
    #[command(verbatim_doc_comment)]
    Flip(FlipArgs),
    /// Compare two BAMs derived from the same reads in the same order, reporting
    /// per-mapQ concordance and per-read discordances.
    ///
    /// Output is TAB-delimited; the first column is a one-letter line type:
    ///
    ///   Q  <mapQ> <a_reads> <a_diff> <a_unmap> <b_reads> <b_diff> <b_unmap> <reads> <diff> <unmapped>
    ///        per-read placement summary. The A group is binned by A's mapQ:
    ///        a_reads = reads mapped in A at this mapQ, a_diff = of those, mapped
    ///        in B but with reciprocal overlap of the exon-block sets < --min-overlap,
    ///        a_unmap = of those, unmapped in B. The B group is the mirror, binned by
    ///        B's mapQ. The last three are binned by q = max(mapQ_A, mapQ_B).
    ///   I  <mapQ> ... (same 9 columns as Q)
    ///        per-read intron-chain summary over SPLICED reads only (>=1 N junction);
    ///        "diff" here means the two junction chains are not identical (same-contig
    ///        required). The trailing trio counts reads spliced in either A or B.
    ///   J  <mapQ> <a_reads> <a_diff> <a_unmap> <b_reads> <b_diff> <b_unmap>
    ///        per-junction summary: a_reads = junctions in A at this (read) mapQ,
    ///        a_diff = of those, no exactly-matching junction in a mapped B, a_unmap =
    ///        of those, B unmapped. B group mirrors it. No max-binned trio.
    ///   U  <#reads>
    ///        pairs unmapped in both files.
    ///   E  <name> <a_ctg> <a_start> <a_end> <a_strand> <a_mapQ> <b_ctg> <b_start> <b_end> <b_strand> <b_mapQ>
    ///        one per discordant pair (only with -e); coordinates are 1-based
    ///        inclusive, strand is +/-, an unmapped end has "." in all five fields.
    #[command(verbatim_doc_comment)]
    Cmp(CmpArgs),
}

#[derive(Args)]
struct FlipArgs {
    /// Input SAM/BAM; "-" or omitted reads from stdin.
    input: Option<PathBuf>,
    /// Output BAM; "-" or omitted writes to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Args)]
struct CmpArgs {
    a: PathBuf,
    b: PathBuf,
    /// Reciprocal-overlap threshold for concordance.
    #[arg(short = 'l', long, default_value_t = 0.5)]
    min_overlap: f64,
    /// Write output to this file instead of stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Emit per-read E lines for discordant reads (off by default).
    #[arg(short = 'e', long = "emit-e")]
    emit_e: bool,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Flip(args) => flip::run(args.input.as_deref(), args.output.as_deref()),
        Command::Cmp(args) => cmp::run(
            &args.a,
            &args.b,
            args.min_overlap,
            args.output.as_deref(),
            args.emit_e,
        ),
    }
}
