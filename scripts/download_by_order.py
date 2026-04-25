#!/usr/bin/env python3
"""
download_by_order.py

Driver that downloads Gammaproteobacteria in chunks, one file per taxonomic
order. Each chunk is a separate download — if one is interrupted, rerun just
that one; the others aren't affected.

Concatenate the resulting files with:
    cat data/by_order/*.fasta.gz > data/gammaproteobacteria_dark.fasta.gz

(gzip files concatenate directly — MMseqs2 and gunzip handle multi-member
gzip transparently.)

Usage:
    python download_by_order.py
    python download_by_order.py --out-dir data/by_order
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Each entry: (output filename stem, NCBI taxonomy ID, human-readable label)
# These are the major orders within Gammaproteobacteria (class 1236).
SUBSETS = [
    ("enterobacterales",   91347,  "Enterobacterales"),
    ("pseudomonadales",    72274,  "Pseudomonadales"),
    ("vibrionales",       135623,  "Vibrionales"),
    ("alteromonadales",   135622,  "Alteromonadales"),
    ("oceanospirillales", 135619,  "Oceanospirillales"),
    ("legionellales",     118969,  "Legionellales"),
    ("xanthomonadales",   135614,  "Xanthomonadales"),
    ("thiotrichales",     72273,   "Thiotrichales"),
    ("chromatiales",      135613,  "Chromatiales"),
    ("aeromonadales",     135624,  "Aeromonadales"),
    ("pasteurellales",    135625,  "Pasteurellales"),
]

# Possible expansion to include all Proteobacteria


# Path to your per-query download script. This driver shells out to it once
# per subset, overriding TAXONOMY_IDS via a small wrapper below.
DOWNLOAD_SCRIPT = Path(__file__).parent / "batch_download_2.py"


def run_one(tax_id: int, label: str, output_path: Path) -> bool:
    """Run the download script for one taxonomic order.

    Returns True on success, False on failure — so the driver can continue
    to the next subset even if one fails.
    """
    if output_path.exists():
        print(f"  SKIP  {label}: {output_path} already exists")
        return True

    print(f"\n=== {label} (tax_id={tax_id}) ===")

    env_override = {"UNIPROT_TAXONOMY_ID": str(tax_id)}

    import os
    env = {**os.environ, **env_override}

    try:
        subprocess.run(
            [sys.executable, str(DOWNLOAD_SCRIPT), "-o", str(output_path)],
            env=env,
            check=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  FAILED {label}: exit code {e.returncode}", file=sys.stderr)
        # Remove partial file so a rerun starts fresh for this subset.
        if output_path.exists():
            output_path.unlink()
        return False
    except KeyboardInterrupt:
        print(f"\n  INTERRUPTED during {label}", file=sys.stderr)
        if output_path.exists():
            output_path.unlink()
        raise


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download Gammaproteobacteria split by taxonomic order.",
    )
    parser.add_argument(
        "--out-dir",
        default="data/by_order",
        help="Directory to write per-order FASTA files (default: %(default)s)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    successes: list[str] = []
    failures: list[str] = []

    for stem, tax_id, label in SUBSETS:
        output_path = out_dir / f"{stem}.fasta.gz"
        if run_one(tax_id, label, output_path):
            successes.append(label)
        else:
            failures.append(label)

    print("\n=== summary ===")
    print(f"  succeeded: {len(successes)}")
    for label in successes:
        print(f"    - {label}")
    if failures:
        print(f"  failed: {len(failures)}")
        for label in failures:
            print(f"    - {label}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)