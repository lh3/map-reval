> [!Warning]
> This project is vibe coded with Claude Code.

# map-reval

`map-reval` compares two alignments of the **same read set against the same
reference** and reports how consistently reads are placed, stratified by
mapping quality.

Another key use is checking a single aligner's **strand symmetry**: align the
reads to a reference and to its per-contig reverse complement, flip the latter
back onto the forward strand, and compare. A strand-symmetric aligner agrees;
the residual disagreement is the signal.

## Build

```sh
cargo build --release      # → target/release/map-reval
```

## Subcommands

- **`flip`** — takes an alignment made against the reverse-complemented
  reference and flips every record (position, strand, CIGAR, SEQ/QUAL, MD/MC/SA
  tags, mate fields) back onto the original strand, emitting a BAM directly
  comparable to a forward-reference alignment. Needs no reference FASTA and is a
  true involution.
- **`cmp`** — compares two BAMs of the same reads (in the same order) and prints
  a TAB-delimited, per-mapQ concordance report.

## Strand-symmetry workflow

```sh
# 1. build the per-contig reverse complement of the reference (same names/lengths)
#    e.g. seqtk seq -r ref.fa > ref-rc.fa
# 2. align the same reads to both (any aligner: bwa, minimap2, STAR, ...)
# 3. flip the RC alignment onto forward coordinates, then compare
map-reval flip rc.bam -o rc.flip.bam
map-reval cmp fwd.bam rc.flip.bam
```

`cmp` also works for any two comparable alignments (two aligners, two parameter
settings, …) of the same reads against the same reference:

```sh
map-reval cmp aln1.bam aln2.bam
```

## Output

`cmp` emits line-typed TAB-delimited rows (see `map-reval cmp --help` for the
full column spec):

- `Q` — per-read placement concordance (reciprocal overlap), by mapQ
- `I` — intron-chain concordance for spliced reads
- `J` — per-junction concordance (exact / shifted / gone / unmapped)
- `U` — reads unmapped in both
- `E` / `F` — per-read / per-junction discordance detail (with `-e`)

## Requirements

Both inputs must list the **same reads in the same order** (name-collated, with
unmapped reads emitted, e.g. `samtools collate`). `cmp` scores each read's
**primary** alignment, but a primary is counted concordant if it matches any
alignment (primary or supplementary) of that read in the other file;
secondary alignments are ignored.
