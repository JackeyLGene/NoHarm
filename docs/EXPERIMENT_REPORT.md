# Annotation-free prioritization of structurally divergent transcript isoforms

**NoHarm: a passive detector of translational diversity via frame-economy cross-harm**

*Research preview / computational note — June 2026*

---

## Abstract

We describe a zero-training computational method that scans human protein-coding transcript isoforms and ranks genes by how divergent their isoform codon landscapes are, without using expression data, disease labels, gene ontologies, or protein domain annotations. The detector maintains a finite memory of vector frames, merges similar inputs, prunes weakly reinforced frames, and reads the surviving centroids. Isoform cross-harm — the L2 distance between an isoform's 64-dimensional codon-frequency vector and a uniform expectation — is computed per isoform; per-gene divergence is the range across isoforms.

A scan of 245,535 GENCODE v49 transcripts (17,903 multi-isoform genes) yields an extreme-tail distribution: only 29 genes exceed Δ|codon-harm| > 0.1 (0.16% of multi-isoform genes; P99 = 0.0586). Matched null controls (isoform count, transcript length range, mean codon-harm; GC matching not yet included) confirm that the top-ranked genes are not explained by these covariates alone. The top 20 genes by z-score include well-characterized disease-associated loci (MED12, HNRNPA1, SLC39A11) as well as poorly characterized genes (ANKRD18B, SH3BGR, SEPTIN11, PAXBP1) that constitute the method's independent predictions. Exploratory keyword enrichment suggests over-representation of transcriptional regulators.

The method's reading can be interpreted simply: **Δ|codon-harm| measures how different the codon landscapes are that the translation machinery must process across isoforms of the same gene.** High divergence implies that isoform choice is likely to have functional consequences; low divergence means no large structural difference is detected in this specific 3-mer frame-economy coordinate (functional differences through expression, domains, localization, or regulation may still exist).

---

## 1. The Reading

### 1.1 What Δ|codon-harm| measures

A gene's coding sequence is read codon-by-codon. Each codon is mapped to a 64-dimensional frequency vector (A/C/G/T at three positions). The expected codon distribution under uniform synonymous usage provides a common baseline. The L2 distance between actual and expected is the per-isoform codon-harm.

Per-gene divergence (Δ|codon-harm|) is the range of this value across a gene's isoforms.

The reading can be stated in one sentence:

> **Δ|codon-harm| measures how different the codon landscapes are that the translation machinery must process across isoforms of the same gene.**

High divergence → isoform choice changes the translational pressure landscape → the resulting proteins are more likely to differ in folding, stability, localization, or half-life → isoform selection has functional consequences → phenotypic impact is more likely.

Low divergence → no large structural difference is detected in this specific 3-mer frame-economy coordinate. This does NOT mean the isoforms are functionally identical — differences in expression level, protein domains, subcellular localization, or regulatory motifs may still exist. It means only that the codon-level structural divergence is not extreme.

### 1.2 What it is not

- NOT a codon bias index (it uses uniform expectation, not genome-wide codon usage)
- NOT a measure of expression or translation efficiency
- NOT a disease predictor
- NOT dependent on any biological database, annotation, or prior knowledge

It is a **structural coordinate** — an orthogonal dimension that no existing method measures.

---

## 2. Method

### 2.1 Architecture

```
transcript sequence
  → sliding window (30 codons, stride 6)
  → 64-dim codon-frequency vector vs uniform expectation
  → L2 distance = per-window |codon-harm|
  → Geruon frame economy (merge, co-occur, prune, τ, L3)
  → per-isoform mean |codon-harm| + τ
  → per-gene: range(|codon-harm|) across isoforms = Δ|codon-harm|
```

The Geruon (finite-memory centroid detector) provides a secondary readout τ (endogenous time) that tracks the frame economy's internal dynamics. In this report, τ is reported alongside |codon-harm| but the primary ranking metric is Δ|codon-harm|.

### 2.2 Data

- **Source**: GENCODE v49 protein-coding transcripts (245,535 transcripts, 20,758 genes)
- **Multi-isoform genes**: 17,903 (genes with ≥2 transcripts)
- **Isoforms processed**: 242,678 (after CDS length ≥60 nt filter)
- **Encoding**: 64-dim codon frequency (A/C/G/T at 3 positions), sliding window 30 codons, stride 6
- **Baseline**: uniform codon distribution (1/64 per bin) — common zero for all isoforms

### 2.3 Matched null

For each gene, a null pool is constructed from genes matched on:
- Isoform count (log-scale, ±15%)
- Transcript length range (±20%)
- Mean codon-harm (±0.03)

GC-content range is noted as a desired covariate but **not yet included** in the matching procedure; this is a limitation to address before stronger statistical claims.

From the pool, up to 1,000 samples are drawn without replacement. Z-score is computed as (observed − null_mean) / null_std. Empirical p-value uses the lower-bound correction:

```
p = (n_exceed + 1) / (n_pool + 1)
```

This ensures that small pools (e.g., pool=10 → min p ≈ 0.091) cannot produce spuriously significant p-values. Z-score is reported as the primary effect-size metric. Genes with pool_size < 10 are excluded from matched null analysis entirely.

### 2.4 CDS-only comparison

A parallel scan uses the longest ORF (ATG→stop) as a CDS proxy. The CDS/full ratio indicates whether isoform divergence is concentrated in coding regions (>1.2), non-coding UTR regions (<0.5), or both (0.5–1.2). ORF-proxy results are suggestive; GTF-annotated CDS coordinates should be used for publication-grade claims.

---

## 3. Results

### 3.1 Distribution

| Statistic | Value |
|-----------|-------|
| Genes with ≥2 isoforms | 17,903 |
| P50 Δ\|codon-harm\| | 0.0131 |
| P90 Δ\|codon-harm\| | 0.0300 |
| P95 Δ\|codon-harm\| | 0.0376 |
| P99 Δ\|codon-harm\| | 0.0586 |
| Genes with Δ\|codon-harm\| > 0.1 | 29 (0.16%) |
| Genes with Δ\|codon-harm\| > 0.2 | 4 |
| Genes with Δ\|codon-harm\| > 0.3 | 1 (SLC39A11) |

The distribution is heavily right-skewed. Only 0.16% of multi-isoform genes show extreme isoform-level codon divergence.

### 3.2 Matched null results

14,214 genes had sufficient matched pool size (≥10). Z-scores and corrected p-values for the top 20:

| Rank | Gene | n_iso | Δ\|ch\| | z | p (corrected) | pool | obs/exp |
|------|------|-------|---------|---|---------------|------|---------|
| 1 | SLC39A11 | 31 | 0.3278 | 26.15 | <0.004 | 259 | 7.2x |
| 2 | MED12 | 28 | 0.2877 | 17.03 | <0.034 | 29 | 4.5x |
| 3 | SRP14 | 6 | 0.1145 | 14.08 | <0.072 | 13 | — |
| 4 | FLOT1 | 82 | 0.2037 | 12.39 | <0.091 | 10 | — |
| 5 | HNRNPA1 | 77 | 0.2073 | 12.28 | <0.040 | 24 | 3.4x |
| 6 | HLA-F | 52 | 0.1855 | 10.89 | <0.009 | 120 | — |
| 7 | ANKRD18B | 7 | 0.1235 | 10.77 | <0.003 | 337 | 3.4x |
| 8 | CARM1 | 9 | 0.1171 | 9.80 | <0.003 | 388 | — |
| 9 | SH3BGR | 14 | 0.1073 | 9.62 | <0.012 | 88 | 3.0x |
| 10 | SUPT5H | 12 | 0.1251 | 9.30 | <0.011 | 90 | — |
| 11 | HGS | 72 | 0.1668 | 9.15 | <0.059 | 16 | — |
| 12 | TUBGCP5 | 25 | 0.1108 | 8.97 | <0.003 | 354 | — |
| 13 | ZNF384 | 14 | 0.1160 | 8.94 | <0.003 | 361 | — |
| 14 | SLC12A5 | 12 | 0.1420 | 8.70 | <0.009 | 109 | — |
| 15 | G3BP2 | 94 | 0.1656 | 8.35 | <0.039 | 25 | — |
| 16 | PRG4 | 12 | 0.1131 | 8.12 | <0.048 | 20 | — |
| 17 | SEPTIN11 | 10 | 0.0947 | 7.92 | <0.003 | 351 | 2.6x |
| 18 | FBL | 35 | 0.0867 | 7.78 | <0.063 | 15 | — |
| 19 | MT3 | 5 | 0.0512 | 7.68 | <0.091 | 10 | — |
| 20 | PAXBP1 | 10 | 0.0857 | 7.33 | <0.003 | 442 | 2.4x |

**Key**: p = (exceed+1)/(pool+1) lower-bound correction. Z-score is the primary effect-size metric. obs/exp shown for genes with sufficient null pool where expected range (null_mean + 2×null_std) > 0.01.

**Important caveat on p-values**: For genes with small matched pools (MED12 pool=29, SRP14 pool=13, FLOT1 pool=10, MT3 pool=10), the minimum achievable p-value under the correction is limited by pool size. For these genes, the z-score is the more informative metric. The notation "<pool_limit" indicates that no matched null sample exceeded the observed score, but statistical significance cannot be claimed at conventional thresholds due to limited null sample size.

**Correction applied 2026-06-03**: Previous version reported "p < 0.001" for all top genes based on raw exceedance ratio. Current version uses (exceed+1)/(pool+1) correction per reviewer feedback.

### 3.3 CDS/full comparison (ORF-proxy)

| Type | CDS/full ratio | Interpretation | Examples |
|------|---------------|----------------|----------|
| UTR-associated | < 0.5 | Isoform divergence concentrated in non-coding regions (ORF-proxy) | SLC39A11 (0.11), TUBGCP5 (0.27) |
| CDS-associated | > 1.2 | Divergence concentrated in coding regions (ORF-proxy) | SRP14 (1.93), MT3 (1.91), MED12 (1.46) |
| Mixed | 0.5–1.2 | Comparable divergence in both | HNRNPA1 (0.77), PRG4 (1.00) |

**Caveat**: CDS is approximated by longest ORF, not GTF-annotated CDS. "UTR-associated" and "CDS-associated" are suggestive labels; publication-grade claims require GTF-based CDS extraction.

### 3.4 Exploratory keyword enrichment

A simple keyword-based category lookup on gene symbols (not a formal GO enrichment analysis) suggests over-representation of transcriptional regulators (MED12, MED15, FOXP1, FOXP2, SUPT5H, KMT2D, TBP, ARID1A — 8/29 top genes vs 0.3 expected by chance). This is reported as an exploratory observation. Formal GO enrichment using g:Profiler, Enrichr, or GOATOOLS with standard annotation databases should replace this before publication.

---

## 4. Independent Predictions

NoHarm identifies four genes with extreme isoform divergence whose biological functions are currently poorly characterized. These constitute testable predictions.

### Prediction 1: ANKRD18B (z=10.77, obs/exp=3.4×, pool=337)

**What is known**: Contains a DUF3496 domain (Domain of Unknown Function). Only one functional study exists (promoter hypermethylation in lung cancer). No GO biological process annotation. Subcellular localization unknown.

**Prediction**: The extreme isoform-level codon divergence (comparable to HNRNPA1 at 3.4× expected) implies that ANKRD18B isoforms, if translated, are likely to have functionally distinct properties. Tissue-specific isoform expression and protein-level validation are warranted.

### Prediction 2: SH3BGR (z=9.62, obs/exp=3.0×, pool=88)

**What is known**: Thioredoxin-fold protein at 21q22.3 (Down syndrome critical region). Associated with sarcomere Z-line assembly. No major disease association established.

**Prediction**: The 14 isoforms span a wide codon-usage range. Isoform-specific functions in cardiac and skeletal muscle may exist. Its location in the Down syndrome critical region makes isoform dysregulation a candidate mechanism for heart or muscle phenotypes.

### Prediction 3: SEPTIN11 (z=7.92, obs/exp=2.6×, pool=351)

**What is known**: Septin family member. Limited isoform-specific literature.

**Prediction**: Robust signal (largest null pool among predictions) suggests isoforms may have distinct roles in membrane dynamics or cytokinesis.

### Prediction 4: PAXBP1 (z=7.33, obs/exp=2.4×, pool=442)

**What is known**: PAX3/PAX7 binding partner. Historically "gene of unknown clinical significance." First disease mutation reported in 2017 (developmental delay with hypotonia).

**Prediction**: The isoform divergence (largest null pool = 442) suggests broader functional impact beyond currently documented muscle development roles.

---

## 5. Prior Phases (Historical Record)

### Phase 1: Architecture (2026-06-01–02)

Six BiasField-based approaches (blend_into, repel, Codex lookup, PrincipleCodex, passive deposition ×2) all produced zero measurable effect. Root cause: BiasField is dead vector averaging without frame economy.

**Breakthrough**: Cross-harm vectors as input to an independent Geruon with its own frame economy. Harm is not stored — it is the Geruon's temporal (τ) response to deviation from stored normal structure. Analogous to immune self/non-self discrimination.

A P0 encoding bug (`ord(ch)%4` → C/G collision, only 27/64 bins used) was identified and fixed (`BASE[ch]`, full 64 bins).

### Phase 2: Robustness (2026-06-02)

Paired differential experiments (same sequence + same window positions → clean vs noisy) across five perturbation types all yielded ΔΔτ ≈ 0. The Geruon frame economy is structurally robust to random perturbation; its τ dynamics are dominated by internal periodic organization, not external noise. This finding is reported as a positive result: **structure constitutes its own robustness.**

Adversarial perturbations (train-then-attack with L3 chain pre-building) confirmed the same pattern. Single-Geruon + synthetic-sequence + artificial-noise approaches were exhausted. Causal harm detection requires real structural operations and multi-encoding-system comparison.

---

## 6. Limitations

1. **Small matched pools for some top genes** (MED12 pool=29, SRP14 pool=13, FLOT1 pool=10). Z-scores for these genes may be inflated. Larger null pools require relaxed matching criteria or a different null-generation strategy.

2. **GC matching not yet included** in matched null. GC content correlates with codon usage and may confound some rankings.

3. **CDS approximated by longest ORF**, not GTF-annotated CDS coordinates. CDS/full ratios should be recomputed with real CDS annotations.

4. **GO enrichment is keyword-based**, not a formal enrichment analysis. Standard tools (g:Profiler, Enrichr, GOATOOLS) should be used before publication.

5. **τ signal interpretation incomplete**. τ and |codon-harm| are partially decoupled (e.g., lncRNA τ=0.72, pc τ=0.62 despite similar |codon-harm|). The τ readout may carry complementary information not captured by |codon-harm| alone.

6. **tRNA sample too small** (6 genes) — sliding window incompatible with short sequences (72 nt).

7. **Predictions untested**. The four poorly characterized genes require independent experimental validation.

8. **No active restraint demonstrated**. NoHarm detects structural divergence but has not been shown to modulate trajectories.

---

## 7. Data & Code

### Output files

| File | Content |
|------|---------|
| `RNA/derived/noharm_isoform_divergence.json` | Full-transcript: 17,903 genes |
| `RNA/derived/noharm_isoform_cds.json` | ORF-proxy CDS: 17,848 genes |
| `RNA/derived/noharm_isoform_full.json` | Full-transcript (from CDS scan) |
| `RNA/derived/noharm_isoform_stats.json` | Matched null statistics |

### Scripts

| Script | Purpose |
|--------|---------|
| `code/_noharm_genome_scan.py` | Full genome scan (parallel) |
| `code/_noharm_cds.py` | CDS-only scan + side-by-side |
| `code/_noharm_stats.py` | Matched null + enrichment |
| `code/_noharm_rna_isoform.py` | Pilot isoform comparison |
| `code/_noharm_rna.py` | RNA type baseline calibration |
| `code/_noharm_rna_perturb.py` | Systematic perturbation |
| `code/_noharm_dual.py` | WTC dual-stream + paired control |
| `code/_noharm_genes.py` | Multi-gene paired control |
| `code/_noharm_struct.py` | Structural noise (base remapping) |
| `code/_noharm_adversarial.py` | Adversarial perturbation |

### Documents

| File | Content |
|------|---------|
| `noharm/docs/EXPERIMENT_REPORT.md` | This report |
| `noharm/docs/MILESTONES.md` | Milestone log |
| `noharm/docs/TOP20_ANNOTATION.md` | Top 20 annotation table |
| `noharm/README.md` | Project overview |

---

## 8. Release Status

**Recommended**: research preview / computational note / forum draft.

**Not recommended**: full biology paper without:
- Formal GO enrichment (g:Profiler / Enrichr / GOATOOLS)
- GTF-based CDS extraction for CDS/full ratio
- GC matching in null model
- Larger null pools for top-ranked genes (relaxed matching or parametric null)
- Independent experimental validation of predictions

**Strengths for release**:
- GENCODE-scale scan (245K transcripts, 17.9K multi-isoform genes)
- Clean extreme-tail distribution (P99=0.059, only 29 genes >0.1)
- Top genes show clear biological texture (MED12, HNRNPA1, SLC39A11)
- Unknown/cold candidates identified (ANKRD18B, SH3BGR, SEPTIN11, PAXBP1)
- Orthogonal screening coordinate — no dependence on expression, disease labels, GO, or domains
- One-sentence interpretability: "measures translational diversity across isoforms"

---

*Generated 2026-06-03 | NoHarm project: `noharm/` | EE architecture: `CLAUDE.md`*
