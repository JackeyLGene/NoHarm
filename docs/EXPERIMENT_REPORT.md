# NoHarm v0.2 GENCODE Report

**Warm-start corrected frame-economy coordinates for transcript isoform and
region screening**

Research preview / computational note, June 2026

Repository: https://github.com/JackeyLGene/NoHarm  
Preprint DOI: https://doi.org/10.5281/zenodo.20518088

---

## Abstract

NoHarm is a zero-training scanner for transcript isoform FASTA files. It groups
transcript isoforms by gene and reports a compact coordinate map:

- `ch_range`: static codon-landscape spread against a shared uniform 64-bin
  baseline.
- `merge_range`: frame-economy response spread, measured as the within-gene
  range of centroid-memory merge rates across isoforms.
- `drift_range`, `churn_range`, `tau_range`: additional detector-response
  traces for audit and region-level spectra.

The v0.2 GENCODE scan adds a warm-start correction before processing each
isoform. This removes a cold-start artifact in v0.1 that had suppressed
`merge_range` and made the response coordinate look artificially sparse.

The corrected scan processed 245,535 GENCODE v49 protein-coding transcript
records and 17,903 multi-isoform genes. The static coordinate remains compact
(`ch_range` P99 = 0.0586; 29 genes exceed 0.1). The warm-start response
coordinate is now continuous at genome scale (`merge_range` P50 = 0.2966,
P99 = 0.6939, max = 0.9046). `ch_range` and `merge_range` are distinct but
moderately coupled (Spearman about 0.62).

The conservative result is that NoHarm provides an annotation-free structural
triage map for transcript isoforms and regions. It does not make clinical,
causal, or mechanistic claims.

---

## 1. Why v0.2 Was Needed

The first public implementation exposed two useful coordinates:

```text
ch_range     static codon-landscape spread
merge_range  detector-response spread
```

An audit found that the v0.1 response coordinate was affected by cold-start.
Every isoform began with an empty centroid memory, so early windows were
disproportionately treated as new frames. This artificially suppressed
`merge_rate` and produced a large zero floor in `merge_range`.

v0.2 fixes this by pre-warming the frame economy with 32 uniform-zero vectors
before each isoform's real windows are processed. This asks a clearer question:

> How does the detector process the isoform stream after its memory has reached
> a saturated baseline state?

`ch_range` is unaffected by this correction because it is computed before
frame-economy processing.

---

## 2. Method

### 2.1 Encoding

For each transcript isoform:

1. Preserve frame positions and normalize bases to A/C/G/T.
2. Read non-overlapping 3-mers in frame 0.
3. Apply a 30-codon sliding window with stride 6.
4. Encode each window as a 64-bin 3-mer frequency vector.
5. Subtract a shared uniform baseline of `1/64` per bin.

### 2.2 Static Coordinate

For each isoform:

```text
mean_ch = mean(L2(window residual))
```

For each gene:

```text
ch_range = max(mean_ch) - min(mean_ch)
```

This is a static codon-landscape spread across isoforms.

### 2.3 Warm-Started Frame-Economy Coordinates

For each isoform:

1. Initialize a finite centroid memory.
2. Process 32 uniform-zero vectors as a warm-start baseline.
3. Process the isoform's residual-window stream.
4. Record merge, drift, churn, tau, and novelty traces.

The public response coordinate is:

```text
merge_rate = merged_windows / total_windows
merge_range = max(merge_rate over isoforms) - min(merge_rate over isoforms)
```

Additional traces include `drift_range`, `churn_range`, `tau_range`, and
`novelty_range`. These are detector-response coordinates. They should not be
described as direct evidence that cellular translation machinery uses the same
mechanism.

---

## 3. Full GENCODE v49 Results

Dataset and runtime:

- GENCODE v49 protein-coding transcript FASTA.
- Total FASTA records: 245,535.
- Scored isoforms: 245,528.
- Multi-isoform genes: 17,903.
- Ranking metric: `merge_range`.
- Runtime: 490.876 seconds on the development machine.

### 3.1 Distribution

Warm-start corrected `merge_range` distribution:

| Coordinate | P50 | P90 | P95 | P99 | Max |
|---|---:|---:|---:|---:|---:|
| `merge_range` | 0.296631 | 0.494180 | 0.550859 | 0.693895 | 0.904641 |

Static coordinate from the matched scan:

| Coordinate | P99 | Genes above 0.1 |
|---|---:|---:|
| `ch_range` | 0.0586 | 29 |

Interpretation:

- `merge_range` is no longer sparse after warm-start.
- The old v0.1 zero floor was an instrument-state artifact.
- `ch_range` remains a compact static residue.
- The response coordinate now highlights broad isoform-processing divergence,
  especially in large and structurally complex genes.

### 3.2 Coordinate Relationship After Correction

The warm-start correction changes the headline from "nearly independent" to
"distinct but coupled":

```text
v0.1 cold-start:
merge_range P50 = 0
83% of genes at zero
Spearman ch_range vs merge_range about 0.21

v0.2 warm-start:
merge_range P50 = 0.2966
no zero floor
Spearman ch_range vs merge_range about 0.62
```

This is a better result, not a weaker one. It means the response coordinate is
not a detached artifact. It shares structure with the static codon landscape,
while still exposing detector-level differences that are not identical to
`ch_range`.

---

## 4. Top Warm-Start Response Genes

Top genes by `merge_range`:

| Rank | Gene | Isoforms | `merge_range` | `mean_ch` | Length Range |
|---:|---|---:|---:|---:|---:|
| 1 | RYR3 | 12 | 0.904641 | 0.208323 | 15468 |
| 2 | GRIN2A | 9 | 0.902708 | 0.206568 | 14606 |
| 3 | TNXB | 23 | 0.900449 | 0.213402 | 13800 |
| 4 | ATRX | 15 | 0.891554 | 0.215850 | 11065 |
| 5 | CREBBP | 7 | 0.889952 | 0.203703 | 10690 |
| 6 | FBXL20 | 10 | 0.887899 | 0.195670 | 10229 |
| 7 | ZEB2 | 42 | 0.884037 | 0.204697 | 9483 |
| 8 | NRXN1 | 29 | 0.882825 | 0.205981 | 9275 |
| 9 | TNC | 30 | 0.877394 | 0.205993 | 8431 |
| 10 | FASN | 13 | 0.877270 | 0.215185 | 8388 |
| 11 | SPTAN1 | 34 | 0.875476 | 0.206022 | 8129 |
| 12 | DMD | 42 | 0.873255 | 0.201632 | 13869 |
| 13 | HS1BP3 | 12 | 0.871560 | 0.206263 | 7650 |
| 14 | PLCB1 | 19 | 0.867965 | 0.204609 | 7223 |
| 15 | GRIN2B | 5 | 0.867806 | 0.200071 | 30531 |
| 16 | SMARCA2 | 61 | 0.865584 | 0.213976 | 6961 |
| 17 | KIF5A | 24 | 0.862851 | 0.214359 | 6691 |
| 18 | SLC12A5 | 12 | 0.855610 | 0.235011 | 6066 |
| 19 | TMEM255B | 20 | 0.854925 | 0.221703 | 6044 |
| 20 | NRDE2 | 7 | 0.847253 | 0.201031 | 13876 |

These genes should be read as high-priority structural candidates, not as
validated biological discoveries. Many are large, multi-isoform, neural,
chromatin, extracellular, or disease-associated genes, which makes them useful
manual-review targets.

---

## 5. Structural Coordinate Map

The current map is:

```text
                         high ch_range
                              |
          static-heavy        |        dual-heavy
          candidates          |        candidates
                              |
low response -----------------+----------------- high response
                              |
          background          |        response-heavy
          / no signal         |        candidates
                              |
                         low ch_range
```

The most conservative public claim is:

> NoHarm exposes a static codon-landscape coordinate and warm-start corrected
> detector-response coordinates. These coordinates are distinct but moderately
> coupled and can be used for first-pass isoform and region triage.

The stronger working hypothesis is:

> Genes or regions extreme in one or more coordinates may deserve biological
> follow-up, especially when they survive matched nulls and comparison with
> standard sequence metrics.

---

## 6. Region Spectra And Case Studies

The same gene-level coordinates can be aggregated into region spectra. This
turns NoHarm from a top-gene list into a first-pass region classifier.

Current exploratory directions:

- structural-production loci such as PRB/KRTAP-like regions;
- immune-diversity loci such as MHC-like regions;
- regulatory-flexibility loci such as transcription-factor or developmental
  clusters;
- maintenance/stability loci used in disease-audit case studies.

Correct current reading:

> NoHarm coordinates provide candidate region spectra that can be tested against
> matched loci and known biological annotations.

---

## 7. AD Audit Boundary

AD-associated gene analysis is useful as an exploratory case study, not as a
clinical result.

Current working interpretation:

- raw disease-set enrichment can be inflated by annotation depth and study
  visibility;
- visibility-matched analysis is required before any disease-level claim;
- tau-related, amyloid-processing, immune, and lysosomal genes may occupy
  different response directions within an AD-associated set;
- no individual gene should be presented as a validated NoHarm disease
  discovery without external evidence and multiple-testing correction.

The AD case is therefore best used to demonstrate careful visibility-aware
evaluation, not to market NoHarm as a disease predictor.

---

## 8. Limitations

1. **Matched nulls are still needed for strong claims.** GC, length, isoform
   count, gene density, and annotation visibility should all be controlled.

2. **Response coordinates depend on detector parameters.** Window size, memory
   cap, merge radius, stride, and warm-start length must be stress-tested.

3. **Biological annotation remains post-hoc.** Named gene categories are
   hypothesis-generating texture, not validation.

4. **No clinical claims.** The tool does not diagnose disease, predict disease,
   or establish causality.

5. **No mechanistic translation claim.** The detector response is a structural
   computation over sequence windows, not a direct model of ribosome behavior.

---

## 9. Recommended Public Framing

Good:

> I built a minimal transcript-isoform scanner with a static codon-landscape
> coordinate and warm-start corrected detector-response coordinates. On GENCODE
> v49, these coordinates are distinct but moderately coupled and nominate
> different genes and regions for follow-up. I am looking for validation advice
> and biological criticism.

Avoid:

> This proves translation machinery uses different processing strategies.

Avoid:

> These are disease genes or validated functional discoveries.

---

## 10. Next Controls

For the next public tool release, prioritize:

1. GC-, length-, isoform-count-, and visibility-matched nulls for response
   coordinates.
2. Comparison with CAI, ENC, GC, codon-usage, conservation, and existing
   isoform-prioritization metrics.
3. GTF-based CDS and UTR extraction.
4. Parameter sweeps for `merge_radius`, memory capacity, window size, stride,
   and warm-start length.
5. Region-level matched nulls.
6. Independent biological review of high-response, high-static, and dual-heavy
   candidates.

---

Generated 2026-06-04. This document supersedes the v0.1 cold-start report.
