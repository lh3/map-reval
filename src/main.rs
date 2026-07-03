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
    Flip(FlipArgs),
    /// Compare two BAM/SAM files derived from the same reads (placeholder).
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
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Flip(args) => flip::run(args.input.as_deref(), args.output.as_deref()),
        Command::Cmp(args) => cmp::run(&args.a, &args.b),
    }
}
