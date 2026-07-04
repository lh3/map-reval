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
    ///   Q  <mapQ> <#reads> <#wrong> <#unmapped>
    ///        per-mapQ summary, binned by q = max(mapQ_A, mapQ_B); #wrong counts
    ///        both-mapped pairs whose reciprocal overlap < --min-overlap; #unmapped
    ///        counts pairs with exactly one unmapped end.
    ///   U  <#reads>
    ///        pairs unmapped in both files.
    ///   E  <name> <a_ctg> <a_start> <a_end> <a_strand> <a_mapQ> <b_ctg> <b_start> <b_end> <b_strand> <b_mapQ>
    ///        one per discordant pair; coordinates are 1-based inclusive, strand is
    ///        +/-, and an unmapped end has "." in all five of its fields.
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
    /// Suppress per-read E lines (emit Q summary only).
    #[arg(long = "no-e")]
    no_e: bool,
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
            !args.no_e,
        ),
    }
}
