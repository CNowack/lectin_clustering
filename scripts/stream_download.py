#!/usr/bin/env python3
"""
download_uniprot.py

Download a UniProt query result set via the /uniprotkb/stream endpoint.

The `stream` endpoint returns all matching entries in a single response — no
pagination, no job polling. Well suited for result sets up to ~10 GB. For
larger sets, fall back to the paginated /search endpoint with Link headers.

Run on any machine with open HTTPS (e.g. your MacBook), then rsync the
resulting .fasta.gz into data/ on the cluster.

Usage:
    python download_uniprot.py
    python download_uniprot.py -o data/proteobacteria.fasta.gz
"""

from __future__ import annotations   # Enables `str | None` syntax on Python 3.9

import sys
import requests
import argparse
from pathlib import Path

# =============================================================================
# QUERY PARAMETERS — mirror check_uniprot_count.py so you can preview then run
# =============================================================================

KEYWORD = None
REVIEWED = False
TAXONOMY_IDS = [1236]
LENGTH_RANGE = [80, 500]
REQUIRE_AFDB = True
PFAM_FAMILIES = None

# =============================================================================
# ENDPOINT + DEFAULT OUTPUT
# =============================================================================

# The stream endpoint returns ALL matching entries in one response. We ask for
# gzip compression to keep the download small; `.fasta.gz` can be read directly
# by MMseqs2 and other tools, or gunzipped after transfer.
UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

DEFAULT_OUTPUT = "data/uniprot.fasta.gz"

# =============================================================================
# QUERY BUILDER
# =============================================================================

def build_query(
    keyword,
    reviewed,
    taxonomy_ids,
    length_range,
    require_afdb,
    pfam_families,
) -> str:
    """Assemble a UniProt query string from the parameter variables."""
    clauses = []

    if keyword:
        clauses.append(f"keyword:{keyword}")

    if reviewed is True:
        clauses.append("reviewed:true")
    elif reviewed is False:
        clauses.append("reviewed:false")

    if taxonomy_ids is not None:
        if isinstance(taxonomy_ids, int):
            taxonomy_ids = [taxonomy_ids]

        if len(taxonomy_ids) == 1:
            clauses.append(f"taxonomy_id:{taxonomy_ids[0]}")
        else:
            tax_clause = " OR ".join(f"taxonomy_id:{t}" for t in taxonomy_ids)
            clauses.append(f"({tax_clause})")

    if length_range:
        lo, hi = length_range
        clauses.append(f"length:[{lo} TO {hi}]")

    if require_afdb:
        clauses.append("database:alphafolddb")

    if pfam_families:
        pfam_clause = " OR ".join(f"xref:pfam-{p}" for p in pfam_families)
        clauses.append(f"({pfam_clause})")

    return " AND ".join(clauses)


# =============================================================================
# DOWNLOAD via the /stream endpoint
# =============================================================================

def download_stream(query: str, output_path: str) -> None:
    """Stream the query results to disk in 1 MB chunks.

    `stream=True` prevents `requests` from loading the whole response into
    memory — essential for multi-GB downloads. Memory usage stays flat at
    ~1 MB regardless of total size.

    `format=fasta` requests FASTA output; `compressed=true` asks the server
    to gzip the stream on the wire — typically 3–5x smaller than raw FASTA.
    """
    response = requests.get(
        UNIPROT_STREAM_URL,
        params={
            "query": query,
            "format": "fasta",
            "compressed": "true",
        },
        stream=True,
    )
    response.raise_for_status()     # Raise if UniProt returned 4xx/5xx

    total_bytes = 0
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1 << 20):   # 1 MiB
            if not chunk:
                continue
            f.write(chunk)
            total_bytes += len(chunk)
            # `\r` rewrites the same line so the progress doesn't scroll.
            print(f"  downloaded: {total_bytes / 1e6:.1f} MB", end="\r")

    # Final newline so subsequent prints don't land on the progress line.
    print(f"\nSaved {output_path} ({total_bytes / 1e6:.1f} MB)")


# =============================================================================
# MAIN
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download a UniProt query result set as gzipped FASTA.",
    )
    parser.add_argument(
        "-o", "--output",
        default=DEFAULT_OUTPUT,
        help="Output path for the .fasta.gz file (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Auto-create parent directory if it doesn't exist. `mkdir -p` semantics.
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    query = build_query(
        keyword=KEYWORD,
        reviewed=REVIEWED,
        taxonomy_ids=TAXONOMY_IDS,
        length_range=LENGTH_RANGE,
        require_afdb=REQUIRE_AFDB,
        pfam_families=PFAM_FAMILIES,
    )

    print(f"Query: {query}")
    print(f"Output: {args.output}\n")

    download_stream(query, args.output)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)