use std::path::Path;

use anyhow::bail;

/// Placeholder for the `cmp` subcommand: compare two alignments of the same reads.
pub fn run(_a: &Path, _b: &Path) -> anyhow::Result<()> {
    bail!("cmp: not yet implemented");
}
