#!/usr/bin/env python3
"""
check_uniprot_count.py

Preview how many sequences a UniProt query would return — without actually
downloading any sequence data. Useful for sizing queries before committing
to a multi-GB async download.

Run on any machine with open HTTPS (e.g. your MacBook).

Usage:
    python check_uniprot_count.py
"""

import requests  # For HTTP requests; handles URL encoding automatically

# =============================================================================
# QUERY PARAMETERS — edit these to customize what you're counting
# =============================================================================

# The UniProt keyword accession. KW-0430 = "Lectin" (binds carbohydrates).
# Set to None to omit the keyword filter.
KEYWORD = None

# Whether to restrict to reviewed (Swiss-Prot) entries only.
# True  = curated entries only (small, high-quality set)
# False = include unreviewed TrEMBL entries (the "dark pool")
# None  = no filter on review status
REVIEWED = False

# NCBI taxonomy ID(s) to restrict the search. Options:
#   None        = no taxonomy filter
#   2           = single clade (Bacteria)
#   [2, 2157]   = multiple clades OR'd together (Bacteria + Archaea)
# Common starting points:
#   2       = Bacteria
#   1224    = Proteobacteria
#   1236    = Gammaproteobacteria
#   1239    = Bacillota/Firmicutes
#   286     = Pseudomonas
#   561     = Escherichia
#   662     = Vibrio
#   11906   = Burkholderia
#   445     = Legionella
#   1301    = Streptococcus
#   1279    = Staphylococcus
TAXONOMY_IDS = [91347]

# Sequence length bracket [min, max] in residues. Set to None to skip.
# Lectin domains typically fall in [80, 800]; tighten to [80, 400] for
# single-domain architectures only.
LENGTH_RANGE = [100, 600]

# Require an AlphaFold structural model to exist for the entry.
# True  = only entries with models in AFDB (what Durairaj et al. need)
# False = no structural filter
REQUIRE_AFDB = True

# Pfam family accessions to include (OR'd together). Set to None or [] to skip.
# These are common carbohydrate-binding domain families:
#   PF00139 = legume lectin, PF00337 = galectin,  PF00059 = C-type lectin,
#   PF00652 = ricin B,       PF01419 = jacalin,   PF07367 = PA14.
PFAM_FAMILIES = None  # e.g. ["PF00139", "PF00337", "PF00059"]

# UniProt REST endpoint for searching UniProtKB. We hit this with HEAD requests
# so the server tells us the total match count without sending any sequences.
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


# =============================================================================
# QUERY BUILDER — assembles the parameters above into a UniProt query string
# =============================================================================

def build_query(
    keyword: str | None,
    reviewed: bool | None,
    taxonomy_ids: list[int] | None,
    length_range: list[int] | None,
    require_afdb: bool,
    pfam_families: list[str] | None,
) -> str:
    """Combine the parameter variables into a single UniProt query string.

    UniProt uses a Lucene-like syntax: `field:value` clauses joined with
    AND / OR / NOT. We build a list of clauses and AND them together.
    """
    clauses = []

    if keyword:
        # `keyword:` matches the controlled-vocabulary keyword accession.
        clauses.append(f"keyword:{keyword}")

    if reviewed is True:
        clauses.append("reviewed:true")
    elif reviewed is False:
        # reviewed:false asks for TrEMBL (unreviewed) entries.
        clauses.append("reviewed:false")
    # If reviewed is None, add no clause — both sets included.

    if taxonomy_ids is not None:
        # Accept either a single int or a list of ints. Normalize to a list
        # so the same code handles both.
        if isinstance(taxonomy_ids, int):
            taxonomy_ids = [taxonomy_ids]

        if len(taxonomy_ids) == 1:
            # Single clade — no parentheses needed.
            clauses.append(f"taxonomy_id:{taxonomy_ids[0]}")
        else:
            # Multiple clades — OR them together, wrap in parens so the OR
            # doesn't spill out and combine with sibling AND clauses.
            tax_clause = " OR ".join(f"taxonomy_id:{t}" for t in taxonomy_ids)
            clauses.append(f"({tax_clause})")

    if length_range:
        lo, hi = length_range
        # UniProt range syntax uses square brackets and TO.
        clauses.append(f"length:[{lo} TO {hi}]")

    if require_afdb:
        # `database:alphafolddb` restricts to entries with AFDB cross-references.
        clauses.append("database:alphafolddb")

    if pfam_families:
        # OR the Pfam xrefs together, then AND the whole group with the rest.
        pfam_clause = " OR ".join(f"xref:pfam-{p}" for p in pfam_families)
        clauses.append(f"({pfam_clause})")

    # Join all clauses with AND.
    return " AND ".join(clauses)


# =============================================================================
# COUNT via the `requests` library
# =============================================================================

def count_via_requests(query: str) -> int | None:
    """Ask UniProt how many entries match the query, without downloading them.

    `requests.head` sends an HTTP HEAD request — headers only, no response
    body. `params={..., "size": 0}` tells UniProt to return zero results per
    page; we only care about the total count, which comes back in the
    `x-total-results` response header.

    `requests` handles URL encoding automatically, so we can pass the raw
    query string with spaces and colons.
    """
    response = requests.head(
        UNIPROT_SEARCH_URL,
        params={"query": query, "size": 0},
    )
    # The count comes back as a string header; convert to int for convenience.
    count_str = response.headers.get("x-total-results")
    return int(count_str) if count_str else None


# =============================================================================
# MAIN — build the query, run the count, print the result
# =============================================================================

def main() -> None:
    # Assemble the query from the parameters at the top of the file.
    query = build_query(
        keyword=KEYWORD,
        reviewed=REVIEWED,
        taxonomy_ids=TAXONOMY_IDS,
        length_range=LENGTH_RANGE,
        require_afdb=REQUIRE_AFDB,
        pfam_families=PFAM_FAMILIES,
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