# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

`map-reval` is a CLI for checking a read aligner's **strand symmetry**: align a read set to a reference and to the per-contig **reverse-complement** of that reference, then compare. If the aligner were perfectly strand-symmetric the two alignments would agree after flipping one back; residual disagreement (e.g. repeat reads placed at a different random hit depending on the reference) is the signal of interest.

Rust binary crate, edition 2024. Built on [`noodles`](https://docs.rs/noodles) for SAM/BAM.

## Commands

```sh
cargo build                 # debug build
cargo test                  # run unit tests (transform helpers live in src/flip.rs)
cargo test flip_pos_chr2    # run a single test by name (substring match)
cargo clippy --all-targets  # lint (kept warning-clean)
cargo fmt
```

## Subcommands

- `map-reval flip [-o OUT] [INPUT]` — reads an alignment made against the RC reference and flips every record back onto the original strand, emitting BAM. `INPUT`/`OUT` default to stdin/stdout (`-` is also accepted).
- `map-reval cmp [-l MIN_OVERLAP] [-o OUT] [-e] [-q MIN_MAPQ] <A> <B>` — compares two BAMs of the same reads (typically a forward-ref BAM vs a `flip`ped RC-ref BAM) and reports per-mapQ concordance. `-l` sets the reciprocal-overlap threshold (default 0.5); `-e` additionally emits per-read `E`/`F` lines (off by default); `-q` (default 5) suppresses those E/F lines when `max(a_mapQ, b_mapQ) < q` (a display filter — the Q/I/J/U aggregates ignore it).

## Architecture

- `src/main.rs` — clap (derive) CLI; `enum Command { Flip, Cmp }` dispatch only.
- `src/flip.rs` — the whole `flip` transform plus pure, unit-tested helpers.
- `src/cmp.rs` — the `cmp` comparison plus the `reciprocal_overlap` helper.

### How `flip` works (the core domain logic)

The RC reference is a **pure per-contig reverse complement**: same contig names and `LN` lengths as the forward reference. Therefore every coordinate transform is derivable from the record plus the header `@SQ` lengths — **no reference FASTA is needed** — and `flip` is a true **involution** (`flip(flip(x)) == x`), which is the primary correctness gate.

Per record (see `transform_record`):
- **POS**: `new_pos = L − pos − span + 2` (1-based), where `span` = reference-consuming CIGAR length and `L` = contig length. Implemented as `flip_pos`.
- **Strand**: toggle FLAG `0x10`; toggle `0x20` when the mate is mapped.
- **CIGAR**: reverse the op vector; **SEQ**: reverse-complement; **QUAL**: reverse.
- **Tags**: `MD:Z` reversed+complemented (`reverse_md`), `SA:Z` per-entry flipped (`transform_sa`), `MC:Z` (mate CIGAR) reversed.
- **Mate fields**: recompute `PNEXT` from `MC` (`new = L_mate − pnext − span(MC) + 2`); negate `TLEN`.
- **Unmapped read** placed at a mapped mate's coordinate: re-place via the mate's `MC`.
- **Header**: preserved as-is (input is typically `SO:unsorted GO:query`, keeping mates adjacent and requiring no re-sort); a `@PG` record is appended.

Tag policy: `flip` transforms `MD`/`MC`/`SA` and treats a fixed allowlist (`KNOWN_TAGS`: `NM/AS/XS/MQ/ms/md/RG/PG/NH/HI` etc.) as orientation-invariant. The **first** record carrying a tag outside that set prints one stderr warning and then passes it through silently (subsequent unknowns are not re-warned).

### How `cmp` works

`cmp` assumes A and B contain the **same reads in the same order** (both must emit unmapped records). It reads one **query-name group** at a time from each file (`next_group`, relying on `GO:query`); a group name mismatch, group-count mismatch, or per-segment presence mismatch is a **fatal error** (no resync). Within a group it buckets records by segment (read1/read2/unpaired), separating each segment's **primary** from its **supplementary** alignments and **ignoring secondaries** (`0x100`). `Placement::of` parses each alignment into a **set of exon blocks** split at `N` plus its **intron junctions** (`ref_blocks`). Same-contig is required for any overlap/match.

**Concordance is supplementary-aware for Q and I.** A read's **primary** interval in A is concordant if it reaches the threshold with *any* alignment (primary or supplementary) of that read in B — and symmetrically (`Counters::compare`). So `a_diff` (A primary vs B's alignment set) and `b_diff` (B primary vs A's set) are independent/directional; a pair is discordant (trio `diff`, E emission) when **either** direction is unmatched. The **J** line stays primary-vs-primary. mapQ binning always uses the primary's mapQ.

Output is TAB-delimited with a one-letter line-type column (documented in full in `map-reval cmp --help`). All summary lines share the `Group { reads, diff, unmap }` shape, one column group binned by A's mapQ, one by B's, and (for Q/I) a trailing trio binned by `max(mapQ_A, mapQ_B)`:
- `Q` — per-read **placement** concordance: `a_diff` = A primary reaches reciprocal overlap `≥ --min-overlap` (default 0.5) with **no** B alignment (primary or supplementary); `unmap` = one end unmapped.
- `I` — per-read **intron-chain** concordance over **spliced reads only** (`a_reads`/`b_reads` require ≥1 junction in the primary; the trio counts reads spliced in either). `a_diff` = A primary's junction chain equals **no** B alignment's chain.
- `J` — per-**junction** concordance (no trio), 8 data cols in order `a_at a_shifted a_gone a_unmap` (+ B mirror): `a_shifted` = no exact match but overlapping a B junction; `a_gone` = B mapped with no overlapping junction; `a_unmap` = B unmapped (exact matches = `a_at − a_shifted − a_gone − a_unmap`).
- `U <#reads>` — pairs unmapped in **both** files.
- `E ...` — one per discordant pair (placement sense), streamed before the summary blocks; read name carries a `/1`/`/2` segment suffix; read extent as a **0-based half-open (BED)** interval, `.` for unmapped ends. Emitted only with `-e`.
- `F ...` — same 12 columns as `E` but per **junction**: one line per non-exact junction on both sides; the focus junction fills its own side (BED interval), the other side shows the largest-overlapping junction (shifted) or `.`×5 (gone / other read unmapped). Emitted only with `-e`. Identical F lines are deduplicated, so a shifted junction (whose A- and B-focused lines coincide) is shown once.

Invariants worth checking: `Σ Q.reads(trio) + U == ` total primary pairs; `Σ I.a_reads` == reads spliced in A; `Σ J.a_reads` == total `N` junctions in A. (`Q.a_diff` and `Q.b_diff` are now directional and need not be equal, since a supplementary can rescue one side but not the other.)

## Testing / verification

Unit tests live next to the code they cover: `src/flip.rs` (`flip_pos`, `revcomp`, `reverse_cigar_str`, `cigar_ref_span`, `reverse_md`, `transform_sa`) and `src/cmp.rs` (`ref_blocks`, `reciprocal_overlap` over block sets). Spliced RNA-seq test BAMs (`HG002.RNA-100k.hs38{f,r}.bam`, minimap2 `splice:sr`) exercise the I/J paths; see [[test-data-bams]].

The strongest end-to-end check for `flip` is the **involution**: `flip` a BAM twice and the result must be byte-identical to the input (`samtools view` diff). For `cmp`, run `flip` on the RC-ref BAM then `cmp` it against the forward-ref BAM and confirm the invariants above; forward-vs-flipped disagreement is expected (highest at low mapQ) — that residual is exactly what the tool measures, not a bug.
