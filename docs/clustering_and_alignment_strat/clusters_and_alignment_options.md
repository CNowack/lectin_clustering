# Handling Overlap Between Lectin Anchors and the Dark Pool

## The Problem

When building the lectin clustering pipeline, the dark pool was downloaded
with `reviewed:false AND taxonomy:Proteobacteria`, while the lectin anchors
were downloaded with the lectin keyword filter. These queries overlap: most
or all of the unreviewed lectins also satisfy the "any unreviewed
Proteobacteria sequence" criterion and ended up in both files.

Concatenating these files naively means lectin sequences appear twice in
the input: once tagged by proteobacterial class (e.g. `dark_gamma`), and
once tagged by review status (e.g. `lectin_rev` / `lectin_unrev`). With
header tagging, a single sequence can only carry one tag, so the question
becomes: which tag matters more for the analysis, and how do we ensure
lectin tags survive the clustering step?

The downstream goal is to find communities with high lectin density —
"hits" being clusters that contain many lectin anchors and unannotated
sequences together, suggesting those unannotated sequences may also be
lectins. For this analysis to work, lectin tags must reliably survive
into the network, and lectin counts per community must be interpretable.

---

## How MMseqs2 Picks Representatives

Understanding representative selection is essential to evaluating each
option, because clustering is where lectin tags can be silently lost.

`mmseqs easy-cluster` runs a multi-step pipeline: a fast k-mer prefilter
to find candidate similar pairs, a Smith-Waterman-style alignment to score
them, then a clustering step that assigns sequences to groups. The
clustering step is where representatives are picked.

The default mode (`--cluster-mode 0`, "set cover") works conceptually
like this: MMseqs2 ranks sequences by their connectivity in the alignment
graph — how many other sequences they have valid alignments to, weighted
by alignment quality. The most-connected sequence becomes a cluster
representative, and every sequence that aligns to it (above the identity
and coverage thresholds) joins that cluster. The next most-connected
sequence among the remaining ones becomes the next representative, and
so on, until every sequence is assigned.

There is a tiebreaker layer: when multiple sequences have similar
connectivity, MMseqs2 prefers longer sequences as representatives. This
is deliberate — longer sequences are more likely to capture full-domain
alignments rather than fragments.

Selection criteria, in order:
1. Sequences with the most connections (highest-degree nodes)
2. Among similarly-connected sequences, the longest one
3. Ties broken by internal ordering

Crucially, the algorithm has no concept of biological importance. It does
not know which sequences are reviewed, well-annotated, or known lectins.
Selection is purely structural — based on graph topology and length.

`--cluster-mode 1` ("connected component") and `--cluster-mode 2`
("greedy") use slightly different rules but share the same core property:
representatives are picked algorithmically, not biologically.

---

## Three Options for Handling the Overlap

### Option 1: Deduplicate, Lectins Only in the Lectin Bucket

Remove every lectin accession from the dark file before clustering. Each
sequence appears exactly once, tagged either as a lectin (`lectin_rev` /
`lectin_unrev`) or as a dark sequence (by proteobacterial class).

**Reasoning to choose this.** Cleanest representation. Each sequence is
one node. When asking "how many lectins are in cluster X?" the count is
unambiguous. Statistical weighting actually means something — a community
with 50 sequences and 3 lectin tags is genuinely 6% lectin.

**The risk.** Lectins enter the clustering step alongside dark sequences.
Whether a given lectin survives as a representative depends on the
structural criteria above. Reviewed lectins are often well-characterized
model proteins from organisms like E. coli, with many close homologs in
TrEMBL. If a lectin has 50 unreviewed homologs in the dataset and one
homolog happens to be slightly longer or more connected, the homolog gets
picked as representative — and the reviewed lectin disappears from
downstream analysis.

The fraction of lectins lost this way is hard to predict without running
it. It could be 10% or it could be 60%, depending on dataset composition
and how MMseqs2's tiebreakers fall.

**Mitigation paths exist but add engineering.** Using `mmseqs cluster`
(not the `easy-cluster` wrapper) with custom representative selection
logic, or doing a two-pass clustering where lectins are consolidated
first then merged with dark sequences under constraints. Both add real
work and the mitigation is still imperfect.

### Option 2: Keep Duplicates Intentionally

Lectins appear once in the dark pool (tagged by proteobacterial class)
and again in the lectin pool (tagged by review status). After clustering
the dark pool, the unclustered lectin file is concatenated on, leaving
lectins present twice in the all-vs-all input.

**Reasoning to choose this.** Guarantees lectin tags survive clustering.
Even if a lectin gets absorbed into a cluster on the dark side, its
second copy in the lectin file is appended after clustering with its tag
intact — so it shows up as its own node in the network. Clustering
reduces redundancy on the dark side; lectins are protected by duplication.

**Implication for analysis.** The network has duplicate lectin nodes —
one as a representative (or absorbed) on the dark side, one as itself
on the lectin side. They are connected by very high-similarity edges in
the all-vs-all (essentially identity matches). When counting lectins
per community, the count is double — divide by 2 for a true count, or
interpret the duplicate as evidence that a lectin appears in the dark
community at a particular position.

**Pragmatic but noisy.** Trades a small amount of analytical noise for
a strong guarantee that lectin tags don't disappear during clustering.

### Option 3: Strip Lectins from Dark, Don't Cluster Lectins at All

Remove lectins from the dark pool *and* skip clustering them entirely.
The dark-minus-lectin sequences cluster down to representatives at 50%
identity. The lectins are appended afterward as raw, unclustered
sequences. Each lectin is its own node in the all-vs-all network.

**Reasoning to choose this.** Every lectin should be findable
individually because they are anchors — losing any to redundancy collapse
defeats their purpose as controls. The dark side benefits from
clustering because the specific identity of any one E. coli homolog
doesn't matter, only the cluster's existence. Lectins should not be
clustered because every lectin is informationally distinct as an anchor.

**Implication for clustering math.** With ~5,000 lectins clustering down
to perhaps 200 representatives at 50% identity, naive clustering would
lose ~4,800 anchors. Keeping them all uncompressed protects against this
and gives maximum power to find lectin-rich communities. The cost is
~5,000 extra sequences in the all-vs-all — trivial compared to the
~2–4M post-clustering dark representatives.

**Why this is the right call for novel-lectin discovery.** The goal is
to find communities with high lectin density. For density math to be
meaningful: lectin tags must reliably survive into the network, and
counting must be unambiguous.

Option 3 gives both for free:
- Every lectin is guaranteed to be a node (no risk of clustering loss)
- The denominator (total sequences in a community) is interpretable —
  N representatives plus any lectins that joined
- The numerator (lectin nodes in the community) is exact
- Density math is meaningful
- Engineering is trivial — one extra concatenation step

The only "cost" is not compressing lectins, which is not really a cost.
The bottleneck is the 13M dark sequences, not the small lectin set.

---

## Summary: Which Option Implies Which Scientific Question

The right option depends on what counts as a "hit":

- "This dark sequence is in the same community as *any* lectin" — want
  as many lectin nodes as possible. Option 3.
- "This community has high lectin density relative to its size" — want
  clean per-sequence counting. Option 1, with extra work to preserve
  tags through clustering.
- "Pragmatic compromise robust to clustering artifacts" — Option 2.

For novel-lectin discovery via lectin density in clusters, Option 3 is
the strongest match.

---

## A Subtler Tagging Issue

Lectins are tagged by review status (`lectin_rev` / `lectin_unrev`),
while dark sequences are tagged by proteobacterial class (`dark_gamma`,
etc.). If a reviewed lectin is also in the gammaproteobacteria dark
file, it has two valid tags. With header-based tagging, only one can
be carried. The `lectin_rev` tag is almost certainly the more important
one for this analysis, since identifying lectins in clusters is the
whole point.

---

## Conceptual Takeaway

**Clustering is a tool for redundancy reduction.** Use it on data where
redundancy reduction is wanted (the dark pool) and skip it on data where
every individual sequence is needed as an anchor (the lectins).