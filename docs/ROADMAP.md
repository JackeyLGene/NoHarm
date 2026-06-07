# NoHarm Roadmap

This roadmap describes the public tool direction after the v0.2 warm-start
GENCODE scan. It is not a manuscript and does not make clinical, causal, or
mechanistic claims.

## Core Direction

The next NoHarm target is a robust structural triage layer for transcript
isoforms and genomic regions.

The core object is not a single gene score. It is a coordinate profile:

```text
static codon-landscape residue
warm-start detector merge response
centroid-norm drift
frame-count churn
stability trace
planned SHP dual-axis coordinate
matched-null and visibility context
```

This profile can be used for first-pass triage:

- which genes have unusually divergent isoforms;
- which isoform pairs define the strongest contrast;
- which loci show a coherent region spectrum;
- which known region classes provide useful anchors for unknown genes or
  regions.

## Track A: Warm-Start Coordinate Map

Goal:

> Establish a reproducible coordinate map in which static and detector-response
> readouts are distinct but not interchangeable.

Current evidence:

- v0.2 removes the v0.1 cold-start zero floor.
- `merge_range` P50 is now 0.2966 and P99 is 0.6939.
- `ch_range` P99 remains 0.0586, with 29 genes above 0.1.
- `ch_range` and `merge_range` are moderately correlated at genome scale
  (Spearman about 0.62), not identical.

Next controls:

- matched nulls for response metrics;
- CAI/ENC/GC/codon-usage/conservation comparisons;
- parameter sweeps;
- GTF-based CDS/UTR extraction.

## Track B: Gene And Region Recognition

Goal:

> Use NoHarm coordinates to classify genes and genomic regions by structural
> processing profile.

Candidate region classes:

- structural production: PRB/KRTAP-like regions;
- immune diversity: MHC/HLA-like regions;
- regulatory flexibility: transcription-factor and developmental loci;
- maintenance/stability: lysosomal, autophagy, calcium, and homeostasis loci.

The public repo should describe this as region spectroscopy or structural
triage. Strong claims require matched random-region nulls.

Current status:

- GeneGrammar has validated a calibrated SHP dual-axis readout on human CDS/UTR
  sequence regimes.
- The NoHarm CLI has not yet integrated this coordinate.
- v0.3 should expose the dual-axis readout as an ordinary output column rather
  than as a separate research script.

## Track C: Disease-Set Audit

Goal:

> Use disease-associated gene sets as difficult case studies for
> visibility-aware evaluation.

Current working interpretation:

- raw disease-set enrichment can be inflated by annotation depth and study
  exposure;
- no disease-level claim should be made without visibility-aware matching;
- subgroup structure may be more informative than whole-label enrichment;
- individual-gene results are exploratory unless they survive multiple-testing
  correction and external validation.

Allowed public framing:

> Disease gene sets are exploratory audits showing how NoHarm handles
> high-visibility biological categories.

Avoid:

> NoHarm detects or predicts a disease.

## Track D: Region-Spectrum Controls

Goal:

> Turn region spectra from suggestive plots into controlled comparisons.

Minimal tests before strong claims:

- compare each target region against length/GC/gene-density matched random
  regions;
- report confidence intervals for region-level coordinates;
- test whether known region classes separate under held-out controls;
- compare NoHarm coordinates against standard genomic summaries;
- publish small reproducible case-study packages.

## Track E: Protein-Fold Extension

Goal:

> Test whether the same structural readout can compare linear sequence rhythm
> with 3D contact-map rhythm in folded proteins.

Current status:

- Early internal SHP-Fold tests are promising but not release-grade.
- This is a separate research track, not a v0.3 dependency.
- Public wording should say "under research", not "protein folding predictor".

Minimum controls before public claims:

- k=20 amino-acid calibration;
- contact-threshold sensitivity;
- matched protein-length and secondary-structure controls;
- comparison against known fold classes and simple composition baselines.

## Release Gate: v0.3 `--dual`

The next major public release should integrate the validated SHP dual-axis
coordinate from GeneGrammar. v0.2 is the corrected baseline; v0.3 should be the
first release that combines:

- the warm-start response-metric correction;
- the existing static and response coordinates;
- an optional `--dual` mode based on calibrated chroma/rhythm SHP readouts;
- a small reproducible demo showing how the dual-axis coordinate changes region
  or gene prioritization.

Minimum gate before release:

- PRB/KRTAP, HOX, and MHC/HLA should be packaged as a reproducible demo with
  matched background controls;
- matched background regions should define a stable null;
- the dual-axis output should be exposed as ordinary CLI output, not as a
  separate research script;
- public wording should remain tool-focused and avoid unresolved mechanism
  claims.

## Recommended Development Order

1. Freeze the v0.2 warm-start coordinate map as the corrected baseline.
2. Integrate the GeneGrammar SHP coordinate as `--dual`.
3. Add a small PRB/KRTAP, HOX, MHC/HLA, and matched-background demo.
4. Add matched-null and metric-comparison controls.
5. Publish v0.3 with a small reproducible case-study package.
6. Keep SHP-Fold as a separate research branch until protein-specific
   calibration is complete.

## Claim Boundary

NoHarm is currently a structural triage tool. It can nominate genes, isoform
pairs, and regions for follow-up. It does not validate biological mechanism,
diagnose disease, or replace existing transcript annotation workflows.
