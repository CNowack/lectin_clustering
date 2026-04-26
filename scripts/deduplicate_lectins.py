#!/usr/bin/env python3
"""
deduplicate_lectins.py

Remove lectin sequences from the dark pool by matching UniProt accessions.

The lectin file and the dark file overlap (since the dark download
filter `reviewed:false AND taxonomy:Proteobacteria` includes the
unreviewed lectins). For Option 3 of the deduplication strategy, the
dark pool must contain only non-lectin sequences before clustering, so
that lectins can be added back as their own anchored nodes after
clustering.

Approach:
    1. Read all UniProt accessions from the lectin FASTA
    2. Stream through the dark FASTA, emitting only entries whose
       accession is NOT in the lectin set
    3. Output a deduplicated dark FASTA ready for clustering

Run on the cluster after tagging and concatenating both files.

Usage:
    python deduplicate_lectins.py
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

# =============================================================================
# INPUT / OUTPUT FILES
# =============================================================================

# Tagged & concatenated lectin file (produced by the tagging script).
# Headers look like: >lectin_rev_tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI ...
LECTIN_FASTA = "data//lectins/lectins.fasta.gz"

# Tagged & concatenated dark file. Headers like: >dark_gamma_tr|XXXXXX|... ...
DARK_FASTA = "data/proteobacteria/proteo_query.fasta.gz"

# Output: dark file with lectin accessions removed.
OUTPUT_FASTA = "data/proteo_query_dedup.fasta.gz"


# =============================================================================
# HEADER PARSING
# =============================================================================

def extract_accession(header_line: str) -> str | None:
    """Extract the UniProt accession from a tagged FASTA header line.

    Tagged header format:
        >tag_tr|A0A1V1IGJ0|A0A1V1IGJ0_ECOLI Description...
        >tag_sp|P12345|HEMA_HUMAN Description...

    The accession is the second '|'-separated field, regardless of the
    tag prefix or the tr/sp database marker.

    Returns None if the header doesn't match the expected format —
    safer than crashing on a malformed line.
    """
    # Strip leading '>' and trailing newline, take only the ID part
    # (everything before the first space).
    first_token = header_line[1:].split(None, 1)[0]

    # Split on '|' and take the second field. Expect at least 3 fields.
    parts = first_token.split("|")
    if len(parts) < 3:
        return None
    return parts[1]


# =============================================================================
# I/O HELPERS
# =============================================================================

def open_fasta(path: str, mode: str):
    """Open .fasta or .fasta.gz transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# =============================================================================
# STEP 1 — collect lectin accessions
# =============================================================================

def collect_lectin_accessions(lectin_path: str) -> set[str]:
    """Read the lectin file once, return the set of accessions seen."""
    accessions: set[str] = set()
    skipped = 0

    with open_fasta(lectin_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            acc = extract_accession(line)
            if acc is None:
                skipped += 1
                continue
            accessions.add(acc)

    print(f"  collected {len(accessions):,} unique lectin accessions")
    if skipped:
        print(f"  warning: {skipped} headers had unexpected format")
    return accessions


# =============================================================================
# STEP 2 — stream through dark, emit only non-lectin entries
# =============================================================================

def filter_dark_pool(
    dark_path: str,
    output_path: str,
    lectin_accessions: set[str],
) -> tuple[int, int]:
    """Stream the dark file, dropping entries with accessions in the lectin set.

    Returns (kept, removed) entry counts.

    Streaming approach: track the current entry's "keep" decision based
    on its header. Subsequent sequence lines inherit that decision until
    the next '>' header line resets it. This keeps memory flat — we never
    hold more than the current line.
    """
    kept = 0
    removed = 0
    keep_current = False

    with open_fasta(dark_path, "rt") as fin, \
         open_fasta(output_path, "wt") as fout:
        for line in fin:
            if line.startswith(">"):
                acc = extract_accession(line)
                if acc is None or acc in lectin_accessions:
                    keep_current = False
                    removed += 1
                else:
                    keep_current = True
                    kept += 1
                    fout.write(line)
            else:
                # Sequence line — write only if the current entry is being kept.
                if keep_current:
                    fout.write(line)

    return kept, removed


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    # Sanity check inputs exist before doing real work.
    for path in (LECTIN_FASTA, DARK_FASTA):
        if not Path(path).exists():
            print(f"ERROR: input file not found: {path}", file=sys.stderr)
            sys.exit(1)

    Path(OUTPUT_FASTA).parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading lectin accessions from {LECTIN_FASTA}")
    lectin_accessions = collect_lectin_accessions(LECTIN_FASTA)

    print(f"\nFiltering dark pool {DARK_FASTA} → {OUTPUT_FASTA}")
    kept, removed = filter_dark_pool(DARK_FASTA, OUTPUT_FASTA, lectin_accessions)

    print("\n=== summary ===")
    print(f"  dark entries kept:    {kept:,}")
    print(f"  dark entries removed: {removed:,} (lectin accessions)")
    print(f"  output file:          {OUTPUT_FASTA}")

    output_size_mb = Path(OUTPUT_FASTA).stat().st_size / 1e6
    print(f"  output size on disk:  {output_size_mb:.1f} MB")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)