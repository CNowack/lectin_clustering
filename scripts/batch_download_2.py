#!/usr/bin/env python3
"""
download_uniprot.py

Download a UniProt query result set via the paginated /uniprotkb/search
endpoint, following `Link` response headers to walk every page.

Why /search instead of /stream:
    - /stream has a hard 10M-result cap; /search has no cap.
    - /search is resumable: each page is a separate HTTP request, so a
      dropped connection only loses the current page, not the whole job.

Run on any machine with open HTTPS (e.g. your MacBook), then rsync the
resulting .fasta.gz into data/ on the cluster.

Usage:
    python download_uniprot.py
    python download_uniprot.py -o data/gammaproteobacteria.fasta.gz
"""

from __future__ import annotations   # Enables `str | None` on Python 3.9

import argparse
import gzip
import re
import sys
import time
from pathlib import Path
import requests

# =============================================================================
# QUERY PARAMETERS — mirror check_uniprot_count.py so you can preview then run
# =============================================================================

KEYWORD = None
REVIEWED = False
#TAXONOMY_IDS = [1224]              # Proteobacteria
LENGTH_RANGE = [100, 600]
REQUIRE_AFDB = True
PFAM_FAMILIES = None

DEFAULT_OUTPUT = "uniprot_dark.fasta.gz"

# Use this to let script be driven by `download_by_order.py`
import os
_env_tax = os.environ.get("UNIPROT_TAXONOMY_ID")
TAXONOMY_IDS = [int(_env_tax)] if _env_tax else [1236]

# =============================================================================
# ENDPOINT + PAGINATION CONFIG
# =============================================================================

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# UniProt returns up to 500 entries per page. Larger values are silently
# capped to 500, so this is the natural chunk size to request.
PAGE_SIZE = 500

# Retry configuration for transient failures (timeouts, 5xx, connection drops).
MAX_RETRIES = 5
RETRY_BACKOFF_SECONDS = 5    # Exponential: 5, 10, 20, 40, 80

# Small delay between successful pages to avoid hammering the server.
# UniProt doesn't publish strict rate limits but ~0.2s between pages is polite.
INTER_PAGE_DELAY = 0.2


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
# PAGINATION — parse the Link header to find the next page URL
# =============================================================================

# UniProt's Link header looks like:
#   <https://rest.uniprot.org/uniprotkb/search?cursor=abc123&size=500>; rel="next"
# This regex pulls out the URL between < and > when rel="next" is present.
LINK_NEXT_RE = re.compile(r'<([^>]+)>;\s*rel="next"')


def extract_next_url(link_header: str | None) -> str | None:
    """Return the 'next' URL from a Link header, or None if no more pages."""
    if not link_header:
        return None
    match = LINK_NEXT_RE.search(link_header)
    return match.group(1) if match else None


# =============================================================================
# HTTP with retries
# =============================================================================

def get_with_retries(url: str, params: dict | None = None) -> requests.Response:
    """GET a URL with exponential backoff on transient failures.

    Retries on: connection errors, timeouts, 5xx server errors, and 429
    (rate limited). Does NOT retry on 4xx client errors like 400/403/404 —
    those indicate a bad request and won't succeed on retry.
    """
    last_exception: Exception | None = None

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, params=params, timeout=60)

            # Retry on rate-limit and server errors.
            if response.status_code == 429 or response.status_code >= 500:
                wait = RETRY_BACKOFF_SECONDS * (2 ** attempt)
                print(
                    f"\n  HTTP {response.status_code}; retrying in {wait}s "
                    f"(attempt {attempt + 1}/{MAX_RETRIES})",
                    file=sys.stderr,
                )
                time.sleep(wait)
                continue

            response.raise_for_status()
            return response

        except (requests.ConnectionError, requests.Timeout) as e:
            last_exception = e
            wait = RETRY_BACKOFF_SECONDS * (2 ** attempt)
            print(
                f"\n  {type(e).__name__}; retrying in {wait}s "
                f"(attempt {attempt + 1}/{MAX_RETRIES})",
                file=sys.stderr,
            )
            time.sleep(wait)

    raise RuntimeError(
        f"Failed after {MAX_RETRIES} retries: {last_exception}"
    )


# =============================================================================
# DOWNLOAD — paginated walk, writing gzipped FASTA incrementally
# =============================================================================

def download_paginated(query: str, output_path: str) -> None:
    """Walk every page of results and append to a gzipped FASTA file.

    We open the output file once with gzip.open() and append each page's
    FASTA bytes as we receive them. This keeps memory flat — we never hold
    more than one page of data at a time.

    Progress reporting uses the x-total-results header from the first page
    to show a percentage.
    """
    # Initial request — subsequent requests follow Link headers with their
    # own cursors, so we only pass `query`/`size`/`format` on page 1.
    next_url: str | None = UNIPROT_SEARCH_URL
    initial_params: dict | None = {
        "query": query,
        "format": "fasta",
        "size": PAGE_SIZE,
        # We do NOT pass compressed=true here — the server's gzip stream is
        # per-page, but we want a single gzip file that concatenates all
        # pages. Easier to gzip on our end as we write.
    }

    total_results: int | None = None
    entries_downloaded = 0
    page_number = 0
    start_time = time.time()

    # gzip.open in 'wb' mode: we write raw FASTA bytes and gzip handles
    # compression on the fly. The resulting .fasta.gz is a valid single-member
    # gzip file readable by any standard tool (MMseqs2, gunzip, zcat).
    with gzip.open(output_path, "wb") as out_file:
        while next_url is not None:
            page_number += 1

            # Only the first request carries query params; follow-ups use
            # the full URL from the Link header (which already has the cursor).
            response = get_with_retries(
                next_url,
                params=initial_params if page_number == 1 else None,
            )
            initial_params = None   # Clear after first use

            # On the first page, grab the total count for progress display.
            if total_results is None:
                total_str = response.headers.get("x-total-results")
                total_results = int(total_str) if total_str else 0
                print(f"  total entries to download: {total_results:,}")

            # Write the page's FASTA bytes straight to the gzipped output.
            # response.content is the full body (safe here — one page is
            # small, a few hundred KB at most).
            out_file.write(response.content)

            # Count entries in this page by counting FASTA headers ('>' at
            # start of line). This is exact and cheap.
            page_entries = response.text.count("\n>") + (
                1 if response.text.startswith(">") else 0
            )
            entries_downloaded += page_entries

            # Progress line, overwritten each page with \r.
            elapsed = time.time() - start_time
            rate = entries_downloaded / elapsed if elapsed > 0 else 0
            pct = (
                f"{100 * entries_downloaded / total_results:.1f}%"
                if total_results
                else "?"
            )
            print(
                f"  page {page_number}: "
                f"{entries_downloaded:,}/{total_results:,} ({pct})  "
                f"{rate:.0f} entries/s",
                end="\r",
            )

            # Find the next page URL in the Link header.
            next_url = extract_next_url(response.headers.get("link"))

            # Brief pause between pages to avoid hammering the server.
            if next_url:
                time.sleep(INTER_PAGE_DELAY)

    # Final summary line (overwrite the progress line).
    elapsed = time.time() - start_time
    size_mb = Path(output_path).stat().st_size / 1e6
    print(
        f"\nSaved {output_path} "
        f"({entries_downloaded:,} entries, {size_mb:.1f} MB, "
        f"{elapsed / 60:.1f} min)"
    )


# =============================================================================
# MAIN
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download a UniProt query result set as gzipped FASTA "
                    "via the paginated /search endpoint.",
    )
    parser.add_argument(
        "-o", "--output",
        default=DEFAULT_OUTPUT,
        help="Output path for the .fasta.gz file (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

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

    download_paginated(query, args.output)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)