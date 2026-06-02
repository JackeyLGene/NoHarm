# Top 20 Isoform-Divergent Genes

Matched-null ranking and biological annotation for the June 2026 research
preview.

Machine-readable table: `data/gencode_v49_top20.tsv`.

## Ranking

The table below is ranked by matched-null z-score. The empirical p-value column
uses the conservative lower-bound correction:

```text
p = (n_exceed + 1) / (pool_size + 1)
```

When no matched-null sample exceeds the observed score, the table reports the
minimum achievable p-value as `<1/(pool+1)`. For small pools, z-score should be
read as the main effect-size statistic rather than as conventional significance.

| Rank | Gene | n isoforms | Delta \|ch\| | z | corrected p | pool | CDS/full | Category |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| 1 | SLC39A11 | 31 | 0.3278 | 26.15 | <0.004 | 259 | 0.11 | Zinc transporter ZIP11 |
| 2 | MED12 | 28 | 0.2877 | 17.03 | <0.034 | 29 | 1.46 | Mediator complex |
| 3 | SRP14 | 6 | 0.1145 | 14.08 | <0.072 | 13 | 1.93 | Signal recognition particle |
| 4 | FLOT1 | 82 | 0.2037 | 12.39 | <0.091 | 10 | 0.30 | Flotillin / lipid raft |
| 5 | HNRNPA1 | 77 | 0.2073 | 12.28 | <0.040 | 24 | 0.77 | Core splicing factor |
| 6 | HLA-F | 52 | 0.1855 | 10.89 | <0.009 | 120 | 0.38 | MHC class I |
| 7 | ANKRD18B | 7 | 0.1235 | 10.77 | <0.003 | 337 | 1.08 | Poorly characterized |
| 8 | CARM1 | 9 | 0.1171 | 9.80 | <0.003 | 388 | 0.30 | Arginine methyltransferase |
| 9 | SH3BGR | 14 | 0.1073 | 9.62 | <0.012 | 88 | 0.65 | Poorly characterized / muscle |
| 10 | SUPT5H | 12 | 0.1251 | 9.30 | <0.011 | 90 | 1.25 | Transcription elongation |
| 11 | HGS | 72 | 0.1668 | 9.15 | <0.059 | 16 | 0.96 | Endosomal signaling |
| 12 | TUBGCP5 | 25 | 0.1108 | 8.97 | <0.003 | 354 | 0.27 | Microtubule complex |
| 13 | ZNF384 | 14 | 0.1160 | 8.94 | <0.003 | 361 | 1.64 | Zinc-finger transcription factor |
| 14 | SLC12A5 | 12 | 0.1420 | 8.70 | <0.009 | 109 | 1.14 | K-Cl cotransporter |
| 15 | G3BP2 | 94 | 0.1656 | 8.35 | <0.039 | 25 | 1.12 | Stress granule assembly |
| 16 | PRG4 | 12 | 0.1131 | 8.12 | <0.048 | 20 | 1.00 | Lubricin |
| 17 | SEPTIN11 | 10 | 0.0947 | 7.92 | <0.003 | 351 | 0.82 | Poorly characterized septin |
| 18 | FBL | 35 | 0.0867 | 7.78 | <0.063 | 15 | 1.23 | Fibrillarin / rRNA processing |
| 19 | MT3 | 5 | 0.0512 | 7.68 | <0.091 | 10 | 1.91 | Metallothionein 3 |
| 20 | PAXBP1 | 10 | 0.0857 | 7.33 | <0.003 | 442 | 1.08 | PAX-binding protein |

## CDS/full Ratio

The CDS/full ratio is computed using a longest-ORF proxy, not GTF-annotated CDS
coordinates. These labels are suggestive and should be recomputed with real CDS
annotations before stronger biological claims.

| Type | CDS/full | Working interpretation | Examples |
|---|---:|---|---|
| UTR-associated | <0.5 | Divergence appears concentrated outside the ORF proxy | SLC39A11, TUBGCP5, FLOT1, CARM1, HLA-F |
| ORF-associated | >1.2 | Divergence appears concentrated inside the ORF proxy | SRP14, MT3, ZNF384, MED12, SUPT5H |
| Mixed | 0.5-1.2 | Comparable divergence in both coordinates | HNRNPA1, HGS, PRG4, PAXBP1 |

## Biological Texture

The top tail contains known biologically structured loci:

- **SLC39A11 / ZIP11**: zinc transporter; strongest full-transcript divergence;
  ORF-proxy comparison suggests the signal is mostly outside the ORF-proxy
  coordinate.
- **MED12**: mediator complex subunit; ORF-associated divergence; known disease
  and drug-resistance relevance.
- **SRP14**: signal-recognition particle component; ORF-associated divergence.
- **FLOT1**: lipid-raft scaffold protein; large isoform set and UTR-associated
  divergence.
- **HNRNPA1**: core RNA-binding and splicing factor; mixed divergence; useful
  calibration point because it is itself a splicing-regulatory node.

The same top tail also contains under-characterized candidates:

- **ANKRD18B**
- **SH3BGR**
- **SEPTIN11**
- **PAXBP1**

These are not claimed as disease genes. They are candidate genes prioritized by
an annotation-free isoform-divergence coordinate and require independent
experimental validation.

## Caveats

- Matched null currently controls isoform count, transcript length range, and
  mean codon-harm. GC matching is not yet included.
- Small matched pools limit p-value resolution for several top genes.
- Keyword enrichment is exploratory and should be replaced by formal GO
  enrichment before publication-grade claims.
- This table is for hypothesis generation, not clinical interpretation.
