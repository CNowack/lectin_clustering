#!/usr/bin/env python3
"""
tag_and_concat.py

Prepend a category tag to each FASTA header, then concatenate all tagged
files into a single combined FASTA. Tags travel with sequences through
MMseqs2 and downstream steps so we can identify reviewed lectins,
unreviewed lectins, and dark sequences in the alignment output.

Tag format: prepended to the sequence ID so MMseqs2 preserves it.
    Original: >tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI ...
    Tagged:   >dark_gamma_tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI ...

Run this once on your collection of input FASTAs. Output is a single
gzipped FASTA ready for MMseqs2.

Usage:
    python tag_and_concat.py
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

# =============================================================================
# INPUT FILES — map each input to its category tag
# =============================================================================

# Map of input file path -> tag prepended to every FASTA header in that file.
# Tags should be short, distinctive, and contain no spaces or pipes (since
# pipes are used as separators in UniProt accessions).
INPUTS = {
    # Lectin anchors — known/curated to be carbohydrate-binding
    # "data/lectins_rev.fasta":      "lectin_rev",
    # "data/lectins_nr.fasta":    "lectin_nr",

    # Dark pool — by Proteobacteria class
    "data/by_order/aeromonadales.fasta.gz":         "dark_aero",
    "data/by_order/alteromonadales.fasta.gz":       "dark_altero",
    "data/by_order/chromatiales.fasta.gz":          "dark_chroma",
    "data/by_order/enterobacterales.fasta.gz":      "dark_entero",
    "data/by_order/legionellales.fasta.gz":         "dark_legion",
    "data/by_order/oceanospirillales.fasta.gz":     "dark_oceano",
    "data/by_order/pasteurellales.fasta.gz":        "dark_pasteur",
    "data/by_order/pseudomonadales.fasta.gz":       "dark_pseudo",
    "data/by_order/thiotrichales.fasta.gz":         "dark_thio",
    "data/by_order/vibrionales.fasta.gz":           "dark_vibrio",
    "data/by_order/xanthomonadales.fasta.gz":       "dark_xantho",
}

OUTPUT_PATH = "data/input_query.fasta.gz"


# =============================================================================
# CORE LOGIC
# =============================================================================

def open_input(path: str):
    """Open .fasta or .fasta.gz transparently for line-by-line text reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")        # 'rt' = read text, decompresses on the fly
    return open(path, "r")


def tag_and_write(input_path: str, tag: str, out_handle) -> int:
    """Read one FASTA, prepend `tag_` to each header, write to out_handle.

    Returns the number of sequences processed.

    Header rewriting:
        >tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI desc...
    becomes:
        >{tag}_tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI desc...
    """
    n_seqs = 0
    with open_input(input_path) as f:
        for line in f:
            if line.startswith(">"):
                # Insert tag right after the '>'. The '_' separator is a
                # safe choice since UniProt IDs use '|' but never '_' at
                # the start.
                out_handle.write(f">{tag}_{line[1:]}".encode("utf-8"))
                n_seqs += 1
            else:
                out_handle.write(line.encode("utf-8"))
    return n_seqs


def main() -> None:
    out_path = Path(OUTPUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # gzip output stream — one combined file, one gzip member.
    total_seqs = 0
    per_file_counts: dict[str, int] = {}

    with gzip.open(out_path, "wb") as out:
        for input_path, tag in INPUTS.items():
            if not Path(input_path).exists():
                print(f"  SKIP   {input_path}: file not found", file=sys.stderr)
                continue

            print(f"  tagging {input_path}  →  prefix '{tag}_'")
            n = tag_and_write(input_path, tag, out)
            per_file_counts[input_path] = n
            total_seqs += n
            print(f"           {n:,} sequences")

    print(f"\nWrote {total_seqs:,} sequences to {out_path}")
    print(f"  size on disk: {out_path.stat().st_size / 1e6:.1f} MB")

    print("\nPer-file counts:")
    for path, n in per_file_counts.items():
        print(f"  {n:>12,}  {path}")


if __name__ == "__main__":
    main()