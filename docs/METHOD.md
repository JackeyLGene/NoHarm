# Method Notes

NoHarm's first standalone release exposes only the minimal scan that produced
the first full GENCODE isoform-divergence ranking.

## Object Of Measurement

The scanner does not claim to discover "harm" in the moral or clinical sense.
It measures divergence among transcript isoforms within the same gene in a
fixed sequence coordinate. The intended reading is:

> large isoform divergence may mark translation-facing structural choices worth
> biological follow-up.

## Fixed Encoder

For each transcript:

1. Clean sequence to A/C/G/T.
2. Read non-overlapping 3-mers from frame 0.
3. Use a 30-codon sliding window, shrinking the window for short transcripts.
4. Compute a 64-bin frequency vector for each window.
5. Subtract the uniform 64-bin baseline.
6. Use the mean vector norm as the transcript score.

For each gene:

```text
ch_range = max(mean_ch over isoforms) - min(mean_ch over isoforms)
```

The highest and lowest scoring isoforms are reported as the primary contrast
pair. The current public scanner reproduces the primary `ch_range` ranking
coordinate used in the June 2026 GENCODE scan. Full research Geruon tau/L3
dynamics are documented in the EE research code; this standalone repo uses a
lightweight tau trace for reporting.

## What The Score Is Not

- Not sequence identity.
- Not expression abundance.
- Not protein-domain annotation.
- Not conservation.
- Not a disease predictor.
- Not a clinical interpretation.

## Immediate Controls Needed For Publications

Before making strong biological claims, report:

- stratification or matched null by number of isoforms;
- matched null by transcript length range;
- matched null by GC range;
- gene-level rather than pair-count-inflated statistics;
- post-hoc biological annotation clearly separated from discovery.

## Release Documents

- `EXPERIMENT_REPORT.md`: full research-preview report.
- `TOP20_ANNOTATION.md`: corrected top-20 table with p-value caveats.
- `WORKFLOW.md`: one-command workflow and processing diagram.
