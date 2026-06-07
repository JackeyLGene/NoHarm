# Method Notes

NoHarm v0.2 is a deliberately minimal transcript-isoform scanner. It exposes an
structural coordinate map rather than a single score:

- `ch_range`: static codon-landscape residue.
- `merge_range`: primary frame-economy detector response.
- `drift_range`, `churn_range`, `tau_range`, `novelty_range`: additional
  detector-response traces for audits and region-level work.

The public release is intended to be easy to run and easy to criticize. It is
not a clinical tool and not a mechanistic RNA-to-protein translation model.

## Object Of Measurement

NoHarm measures divergence among transcript isoforms within the same gene in a
fixed sequence-derived coordinate system.

The intended reading is:

> Large or unusual isoform divergence may mark genes or loci worth biological
> follow-up.

The current tool does not directly compare gene, transcript, CDS, protein
sequence, and protein structure layers. Later NoHarm work may add multi-layer
consistency checks. A calibrated SHP dual-axis coordinate has been validated in
the companion GeneGrammar work and is the intended basis for the future
`--dual` mode.

### v0.2 Warm-Start Correction

v0.2 adds a mandatory warm-up step: before processing each isoform's real window
vectors, the FrameEconomy processes 32 uniform-zero vectors to fill memory and
reach its saturated processing regime. This eliminates a cold-start artifact
discovered during Phase 5 fork controls. In v0.1, the first windows of every
isoform always encountered an empty or weakly-formed memory, artificially
suppressing `merge_rate` and inflating apparent sparsity (83% of genes had
`merge_range = 0`). The correction brings `merge_range` P50 from 0 to ~0.30.
`ch_range` (static L2) is unaffected by this change as it is computed from
encoder output before the FrameEconomy processes any vectors.

## Fixed Encoder

For each transcript:

1. Preserve frame positions and normalize bases to A/C/G/T.
2. Read non-overlapping 3-mers from frame 0.
3. Use a 30-codon sliding window, shrinking the window for short transcripts.
4. Compute a 64-bin frequency vector for each window.
5. Subtract the uniform 64-bin baseline.

This produces a stream of residual 64-dimensional window vectors.

## Static Coordinate: `ch_range`

For each isoform:

```text
mean_ch = mean(L2(window_vector - uniform_baseline))
```

For each gene:

```text
ch_range = max(mean_ch over isoforms) - min(mean_ch over isoforms)
```

This is the static composition-like coordinate. It is useful as an ablation
baseline and as a candidate generator. It should not be called a codon-bias
metric in the narrow CAI/ENC sense because it does not use synonymous codon
groups, expression references, or the genetic code.

## Frame-Economy Coordinates

Each isoform's window-vector stream is passed through a finite centroid memory:

- memory capacity: 16 frames;
- merge radius: 0.08;
- each incoming vector either merges into the nearest frame or creates a new
  frame;
- if capacity is exceeded, the weakest frame is dropped.

The primary response coordinate is:

```text
merge_rate = merged_windows / total_windows
merge_range = max(merge_rate over isoforms) - min(merge_rate over isoforms)
```

`merge_range` measures how differently the NoHarm detector processes isoforms
from the same gene. It is a detector-response coordinate. It does not prove that
cellular translation machinery itself uses the same strategy.

Additional traces are reported for audits and case studies:

- `tau_range`: range of the stability trace.
- `novelty_range`: range of novelty coefficient-of-variation.
- `drift_range`: range of centroid-norm drift.
- `churn_range`: range of frame-count churn.

These traces can be useful for region spectra and disease-audit studies, but
the public result should not be reduced to any one of them.

## Structural Coordinate Map

The core NoHarm output is a coordinate map:

```text
                         high ch_range
                              |
          static-only         |        dual-extreme
          candidates          |        candidates
                              |
low response -----------------+----------------- high response
                              |
          background          |        response-only
          / no signal         |        candidates
                              |
                         low ch_range
```

In this map:

- static-only candidates differ in codon landscape but do not strongly change
  the detector's processing behavior;
- response-only candidates look quiet to the static coordinate but provoke
  different detector responses across isoforms;
- dual-extreme candidates are high in both coordinates and are the most obvious
  manual-review targets.

This two-coordinate design is the main reason NoHarm should be treated as a
structural screening instrument rather than as a single k-mer statistic.

## Planned Dual-Axis Coordinate

The planned v0.3 extension is a dual-axis SHP coordinate. Unlike `ch_range`,
which summarizes a static 64-bin codon landscape, SHP compares two binary views
of the same local nucleotide window:

- chroma: which 3-mers are present;
- rhythm: which adjacent 3-mer transitions occur;
- cross-harm: Jaccard distance between the two activation sets;
- fixed_wit: event rate above a fair-IID calibration threshold.

The current calibration inherited from GeneGrammar uses:

```text
k = 4 nucleotide alphabet
n = 3 k-mer granularity
D = 64 hash buckets
W = 128 nt windows
theta_0 = 0.0999
```

This coordinate is intended for CDS/UTR, region-spectrum, and gene-regime
screening. It is not yet exposed in the public NoHarm CLI. Until it is
integrated, NoHarm v0.2 should be read as the warm-start corrected static +
frame-economy baseline.

## Region Spectra

Gene-level coordinates can be aggregated across neighboring genes or selected
loci to form a region spectrum. A region spectrum may include:

- mean and spread of `ch_range`;
- mean and spread of response coordinates;
- adjacent-gene similarity;
- outlier genes within a locus;
- centroid-landscape comparisons between loci.

This is useful for exploratory region classification. For example, structural
production loci, immune-diversity loci, and regulatory loci may occupy different
parts of the coordinate map. These are case-study hypotheses, not public
claims.

## Relationship To Existing Metrics

The NoHarm coordinates are distinct from several familiar measures:

- **Not GC content.** The vector spans all 64 trinucleotide bins, not only
  aggregate G+C frequency.
- **Not CAI or ENC.** NoHarm does not consult synonymous codon groups, highly
  expressed reference genes, or the genetic code.
- **Not sequence conservation or dN/dS.** It is computed from transcript
  sequences without alignments or phylogenetic trees.
- **Not expression abundance.** No expression values are used.
- **Not annotation visibility.** No disease labels, GO terms, PubMed counts, or
  protein-domain annotations are used during scanning.

Publication-grade biological claims still require direct comparison with these
standard covariates and metrics.

## Protein Folding Direction

NoHarm currently operates on transcript-derived symbol streams. A separate
protein-folding extension is under research: the same chroma/rhythm idea can be
tested with a protein alphabet and a rhythm axis defined by 3D contact maps
rather than linear sequence adjacency. Early internal SHP-Fold experiments
suggest that the choice of rhythm axis can expose fold-dependent structure that
linear sequence scans miss.

This is a research direction only. It is not part of the v0.2 CLI and should
not be treated as a validated protein-structure predictor.

## Current Empirical Pattern

On the warm-start corrected GENCODE v49 protein-coding transcript FASTA:

- `ch_range` has a compact right tail: P99 = 0.0586, with 29 genes above 0.1.
- `merge_range` is no longer sparse after warm-start: P50 = 0.2966,
  P99 = 0.6939, max = 0.9046.
- No genes remain at the old cold-start zero floor.
- `ch_range` and `merge_range` are distinct but moderately coupled at genome
  scale (Spearman about 0.62 in the current implementation).
- The v0.1 sparse response result (`merge_range` P50 = 0; 83% zero; Spearman
  about 0.21) was inflated by cold-start and should be treated as superseded.

This is the current core result.

## What The Score Is Not

- Not a disease predictor.
- Not a clinical interpretation.
- Not a causal claim.
- Not direct evidence of translation efficiency.
- Not protein-domain annotation.
- Not conservation.
- Not a replacement for existing transcript or variant annotation tools.

## Immediate Controls Needed For Strong Claims

Before making strong biological claims, report:

- GC-, length-, isoform-count-, and visibility-matched nulls;
- matched nulls for response metrics, not only for `ch_range`;
- GTF-based CDS and UTR extraction instead of longest-ORF proxies;
- comparison with CAI, ENC, GC, codon-usage, and existing isoform-prioritization
  metrics;
- gene-level and region-level statistics without pair-count inflation;
- post-hoc biological annotation clearly separated from discovery.

## Release Documents

- `EXPERIMENT_REPORT.md`: current GENCODE report.
- `WORKFLOW.md`: one-command workflow and processing diagram.
- `ROADMAP.md`: current v0.2 direction and case-study boundaries.
- `GENEGRAMMAR_BRIDGE.md`: SHP dual-axis bridge planned for v0.3.
- `TOP20_ANNOTATION.md`: legacy/static-coordinate annotation table.

## Broader Theory

The public tool is intentionally usable as a small bioinformatics scanner.
Readers interested in the frame-economy theory behind the name can browse:
https://jackeylgene.github.io/GBE
