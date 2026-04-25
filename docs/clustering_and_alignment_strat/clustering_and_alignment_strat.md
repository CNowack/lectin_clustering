# Handling Lectin Anchors in the Clustering Pipeline — Chosen Approach

## The Problem in One Sentence

The dark pool (unreviewed Proteobacteria) contains the lectin anchors as
a subset, since the queries overlap. Naive clustering of the combined
set risks losing lectin tags because MMseqs2 picks representatives
algorithmically — based on graph connectivity and length, not on
biological importance.

## How MMseqs2 Picks Representatives (Brief)

`mmseqs easy-cluster` ranks sequences by connectivity in the alignment
graph (most-connected first), with longer sequences preferred as
tiebreakers. The most-connected sequence becomes a representative; every
sequence aligning to it joins that cluster. This continues until all
sequences are assigned. MMseqs2 has no concept of biological importance
— it cannot prefer a reviewed lectin over its unreviewed homologs. If
a homolog happens to be longer or more connected, the homolog wins and
the lectin is absorbed into the cluster, losing its tag from downstream
analysis.

## The Approach

Strip lectin sequences out of the dark pool, then skip clustering the
lectins entirely:

1. Cluster the dark-minus-lectin pool to representatives at 50% identity
   (collapses 13M sequences to perhaps 2–4M)
2. Append all lectins back as raw, unclustered sequences
3. Run all-vs-all on the combined set

## Why This Works for Novel-Lectin Discovery

The downstream goal is to find communities with high lectin density.
For density math to be meaningful, two things must hold: every lectin
must reliably appear as a node in the network, and counts must be
unambiguous.

This approach delivers both:
- Lectins never enter the clustering step, so no risk of tag loss
- Each lectin is guaranteed to be its own node in the all-vs-all
- Community lectin density (lectins / total nodes) is exact and
  interpretable
- Engineering is trivial — one extra concatenation step

The trade-off is not compressing the lectin set. With only ~5,000
lectins, this is negligible compared to the 13M-sequence dark pool that
genuinely needs compression.

## Conceptual Principle

Clustering is a tool for redundancy reduction. Apply it where redundancy
reduction is wanted (the dark pool) and skip it where every individual
sequence is needed as an anchor (the lectins).