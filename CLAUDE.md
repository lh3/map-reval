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
- `map-reval cmp [-l MIN_OVERLAP] [-o OUT] [--no-e] <A> <B>` — compares two BAMs of the same reads (typically a forward-ref BAM vs a `flip`ped RC-ref BAM) and reports per-mapQ concordance. `-l` sets the reciprocal-overlap threshold (default 0.5); `--no-e` suppresses per-read lines.

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

`cmp` assumes A and B contain the **same reads in the same order** (both must emit unmapped records). It iterates **primary alignments only** in lock-step (skipping `0x100`/`0x800`), pairing the k-th primary of each; a read-name/segment mismatch or a primary-count mismatch is a **fatal error** (no resync). Concordance is **reciprocal overlap** of the two reference footprints (`[pos, pos+span)`, same-contig required): `|A∩B| / |A∪B| ≥ --min-overlap` (default 0.5). Splicing (`N`) is folded into the span for now (approximate; warns once). Each read is binned by `q = max(mapQ_A, mapQ_B)`.

Output is TAB-delimited with a one-letter line-type column (documented in full in `map-reval cmp --help`):
- `Q <mapQ> <#reads> <#wrong> <#unmapped>` — per-mapQ summary; `#wrong` = both-mapped but overlap < threshold; `#unmapped` = exactly one end unmapped.
- `U <#reads>` — pairs unmapped in **both** files.
- `E <name> <a_ctg a_start a_end a_strand a_mapQ> <b_...>` — one per discordant pair (streamed before the `Q`/`U` block); unmapped ends use `.` for all five fields. Suppressed by `--no-e`.

Invariant worth checking after changes: `Σ Q.#reads + U == ` total primary pairs, and `#E lines == Σ (Q.#wrong + Q.#unmapped)`.

## Testing / verification

Unit tests live next to the code they cover: `src/flip.rs` (`flip_pos`, `revcomp`, `reverse_cigar_str`, `cigar_ref_span`, `reverse_md`, `transform_sa`) and `src/cmp.rs` (`reciprocal_overlap`).

The strongest end-to-end check for `flip` is the **involution**: `flip` a BAM twice and the result must be byte-identical to the input (`samtools view` diff). For `cmp`, run `flip` on the RC-ref BAM then `cmp` it against the forward-ref BAM and confirm the invariants above; forward-vs-flipped disagreement is expected (highest at low mapQ) — that residual is exactly what the tool measures, not a bug.
