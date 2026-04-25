#!/usr/bin/env python3
"""
check_uniprot_count.py

Preview how many sequences a UniProt query would return — without actually
downloading any sequence data. Useful for sizing queries before committing
to a multi-GB download.

Run on any machine with open HTTPS (e.g. your MacBook).

Usage:
    python check_uniprot_count.py 1224                 # single taxonomy ID
    python check_uniprot_count.py 1224 1236            # multiple OR'd together
    python check_uniprot_count.py 1224 --length 80 800 # override length range
    python check_uniprot_count.py 1224 --no-afdb       # drop AFDB filter
    python check_uniprot_count.py --help               # see all options
"""

from __future__ import annotations    # `str | None` syntax on Python 3.9

import argparse
import requests

# =============================================================================
# DEFAULT QUERY PARAMETERS — overridable via CLI flags
# =============================================================================

# The UniProt keyword accession. KW-0430 = "Lectin" (binds carbohydrates).
DEFAULT_KEYWORD = None

# Review status filter:
#   True  = curated Swiss-Prot only
#   False = unreviewed TrEMBL only (the "dark pool")
#   None  = both
DEFAULT_REVIEWED = False

# Sequence length bracket [min, max] in residues.
DEFAULT_LENGTH_RANGE = [100, 600]

# Require an AlphaFold model to exist for the entry.
DEFAULT_REQUIRE_AFDB = True

# Pfam family accessions (OR'd). Common carbohydrate-binding families:
#   PF00139 = legume lectin, PF00337 = galectin,  PF00059 = C-type lectin,
#   PF00652 = ricin B,       PF01419 = jacalin,   PF07367 = PA14.
DEFAULT_PFAM_FAMILIES = None

# UniProt REST endpoint.
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


# =============================================================================
# COMMON TAXONOMY IDS — quick reference for the CLI
# =============================================================================

TAXONOMY_REFERENCE = """
Common NCBI taxonomy IDs:
  2       Bacteria                    1224    Proteobacteria
  2157    Archaea                     1236    Gammaproteobacteria
  2759    Eukaryota                   1239    Bacillota/Firmicutes
  10239   Viruses                     286     Pseudomonas
  33208   Metazoa (animals)           561     Escherichia
  33090   Viridiplantae (plants)      662     Vibrio
                                      11906   Burkholderia
                                      445     Legionella
                                      1301    Streptococcus
                                      1279    Staphylococcus
"""


# =============================================================================
# QUERY BUILDER
# =============================================================================

def build_query(
    keyword: str | None,
    reviewed: bool | None,
    taxonomy_ids: list[int] | None,
    length_range: list[int] | None,
    require_afdb: bool,
    pfam_families: list[str] | None,
) -> str:
    """Assemble a UniProt query string from the parameter values."""
    clauses = []

    if keyword:
        clauses.append(f"keyword:{keyword}")

    if reviewed is True:
        clauses.append("reviewed:true")
    elif reviewed is False:
        clauses.append("reviewed:false")

    if taxonomy_ids:
        if len(taxonomy_ids) == 1:
            clauses.append(f"taxonomy_id:{taxonomy_ids[0]}")
        else:
            # Parens required so the OR doesn't bleed into sibling ANDs.
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
# COUNT via UniProt's HEAD response header
# =============================================================================

def count_via_requests(query: str) -> int | None:
    """Ask UniProt how many entries match the query, without downloading them.

    HEAD request fetches headers only. `size=0` tells UniProt to return zero
    results per page. The total count comes back in the `x-total-results`
    response header.
    """
    response = requests.head(
        UNIPROT_SEARCH_URL,
        params={"query": query, "size": 0},
    )
    count_str = response.headers.get("x-total-results")
    return int(count_str) if count_str else None


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Preview UniProt query result counts without downloading.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=TAXONOMY_REFERENCE,
    )

    # Positional: one or more taxonomy IDs (OR'd together if multiple)
    parser.add_argument(
        "taxonomy_ids",
        type=int,
        nargs="+",
        help="One or more NCBI taxonomy IDs (multiple are OR'd).",
    )

    # Length range — pass two ints, or pass --no-length to skip
    parser.add_argument(
        "--length",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=DEFAULT_LENGTH_RANGE,
        help=f"Sequence length range in residues "
             f"(default: {DEFAULT_LENGTH_RANGE[0]} {DEFAULT_LENGTH_RANGE[1]})",
    )
    parser.add_argument(
        "--no-length",
        action="store_true",
        help="Disable the length filter entirely.",
    )

    # Review status — mutually exclusive flags map to True/False/None
    review_group = parser.add_mutually_exclusive_group()
    review_group.add_argument(
        "--reviewed",
        action="store_true",
        help="Only reviewed (Swiss-Prot) entries.",
    )
    review_group.add_argument(
        "--unreviewed",
        action="store_true",
        help="Only unreviewed (TrEMBL) entries — default.",
    )
    review_group.add_argument(
        "--all-review",
        action="store_true",
        help="Both reviewed and unreviewed.",
    )

    # AFDB requirement
    parser.add_argument(
        "--no-afdb",
        action="store_true",
        help="Don't require an AlphaFold model (default: require).",
    )

    # Other filters
    parser.add_argument(
        "--keyword",
        type=str,
        default=DEFAULT_KEYWORD,
        help="UniProt keyword accession (e.g. KW-0430 for Lectin).",
    )
    parser.add_argument(
        "--pfam",
        type=str,
        nargs="+",
        default=DEFAULT_PFAM_FAMILIES,
        help="Pfam family accessions to OR together (e.g. PF00139 PF00337).",
    )

    return parser.parse_args()


def resolve_review_status(args: argparse.Namespace) -> bool | None:
    """Map the mutually-exclusive review flags to True/False/None."""
    if args.reviewed:
        return True
    if args.all_review:
        return None
    # Default and explicit --unreviewed both fall here
    return False


def main() -> None:
    args = parse_args()

    # Apply --no-length to override the default range
    length_range = None if args.no_length else args.length

    query = build_query(
        keyword=args.keyword,
        reviewed=resolve_review_status(args),
        taxonomy_ids=args.taxonomy_ids,
        length_range=length_range,
        require_afdb=not args.no_afdb,
        pfam_families=args.pfam,
    )

    print(f"Query: {query}\n")

    count = count_via_requests(query)
    if count is not None:
        # Rough size estimate: ~350 bytes per FASTA record on average.
        estimated_gb = count * 350 / 1e9
        print(f"Total results: {count:,}")
        print(f"Estimated FASTA size: ~{estimated_gb:.2f} GB")
    else:
        print("No count returned — check the query syntax.")


if __name__ == "__main__":
    main()