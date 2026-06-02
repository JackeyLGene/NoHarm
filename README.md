# NoHarm

Annotation-free prioritization of structurally divergent transcript isoforms.

NoHarm is a lightweight bioinformatics scanner for ranking genes whose transcript
isoforms diverge strongly in a fixed 3-mer codon-landscape coordinate. It is
designed as a low-cost triage layer: run it on a transcript FASTA, get a short
list of genes and isoform pairs that may deserve biological follow-up.

The tool does **not** make clinical claims, diagnose disease, or infer causal
mechanisms. It produces a reproducible structural ranking.

## At A Glance

```mermaid
flowchart TD
    A["Transcript isoforms"] --> B["Encode each isoform<br/>as 3-mer/codon windows"]
    B --> C["Compare with a shared baseline<br/>uniform 64-bin expectation"]
    C --> D["Read structural residue<br/>mean |codon-harm| per isoform"]
    D --> E["Compare residues within each gene<br/>max-min divergence"]
    E --> F["Rank candidate genes<br/>and contrast isoform pairs"]
```

In one sentence:

> NoHarm asks how much the codon landscapes of isoforms from the same gene
> differ after being projected onto a shared baseline.

## Why This Exists

Long-read transcriptomics and modern genome annotation produce very large
isoform sets. The hard question is often not "are there isoforms?", but:

> Which genes have isoforms whose codon landscapes diverge enough to deserve a
> closer look?

NoHarm provides one annotation-free coordinate for that triage problem.

## Quick Start

From the repository root:

```bash
python -m pip install -e .
noharm scan --fasta data/demo_isoforms.fa --out results/demo
```

Or without installation:

```bash
$env:PYTHONPATH="src"
python -m noharm scan --fasta data/demo_isoforms.fa --out results/demo
```

For a GENCODE transcript FASTA:

```bash
noharm scan --fasta gencode.v49.pc_transcripts.fa.gz --out results/gencode_v49 --workers 8
```

## Outputs

`noharm scan` writes:

- `gene_divergence.tsv` — gene-level isoform divergence ranking.
- `isoform_scores.tsv` — per-transcript scores.
- `summary.json` — parameters, distribution, and top genes.
- `report.md` — small human-readable report.

Key columns:

- `ch_range`: max isoform mean cross-harm minus min isoform mean cross-harm.
- `ch_std`: within-gene dispersion of isoform scores.
- `tau_range`: spread in the lightweight residue trace.
- `min_ch_transcript`, `max_ch_transcript`: the top contrast pair.

## Research Preview Results

The June 2026 GENCODE v49 scan processed 245,535 protein-coding transcripts,
covering 17,903 multi-isoform genes. The resulting distribution had a small
extreme tail: P99 `Delta |codon-harm| = 0.0586`, with only 29 genes above `0.1`.

Top matched-null genes include known biologically structured loci such as
`SLC39A11`, `MED12`, `SRP14`, `FLOT1`, `HNRNPA1`, and `HLA-F`, plus
under-characterized candidates such as `ANKRD18B`, `SH3BGR`, `SEPTIN11`, and
`PAXBP1`.

See:

- [Experiment Report](docs/EXPERIMENT_REPORT.md)
- [Top 20 Annotation](docs/TOP20_ANNOTATION.md)
- [Top 20 TSV](data/gencode_v49_top20.tsv)
- [Workflow](docs/WORKFLOW.md)
- [Method Notes](docs/METHOD.md)

## Default Method

The public scanner intentionally fixes the encoder:

1. Split each transcript into non-overlapping 3-mer positions.
2. Slide a 30-codon window across the transcript.
3. Convert each window into a normalized 64-dimensional 3-mer vector.
4. Subtract the uniform 64-bin baseline.
5. Average the vector magnitude over windows to score each isoform.
6. Rank each gene by the range of isoform scores.

No GO terms, disease labels, expression values, protein domains, or genetic-code
semantics are used by the default ranking.

## Current Interpretation

NoHarm should be read as:

> An annotation-free prioritization layer for transcript isoform codon-landscape divergence.

The score is a candidate generator. Known high-ranking genes can calibrate the
signal; under-characterized high-ranking genes become follow-up candidates.
Because codon usage can influence downstream biology, high-ranking genes may be
useful candidates for follow-up, but the v0.1 scanner does not directly model
gene-to-protein translation.

## Important Caveats

- The standalone public scanner reproduces the primary `Delta |codon-harm|`
  ranking coordinate. It compares isoforms to a shared uniform baseline; it does
  not yet implement the deeper gene/transcript-to-CDS/protein alignment planned
  for later NoHarm work.
- Matched null controls currently cover isoform count, transcript length range,
  and mean codon-harm. GC matching is not yet included.
- CDS/full comparisons currently use a longest-ORF proxy, not GTF-annotated CDS
  coordinates.
- Biological annotations are for hypothesis generation and require independent
  validation.

## Status

Research preview, June 2026.

Planned next steps:

- Matched null controls by isoform count, transcript length, GC range, and mean score.
- CDS-only and UTR-aware modes.
- Gene-level codex construction from isoform traces.
- Optional protein/topology annotation joins.

## Related Work

NoHarm is a practical tool extracted from a broader frame-economy research
program. Readers interested in the underlying theory can browse the GBE project:
[https://jackeylgene.github.io/GBE](https://jackeylgene.github.io/GBE).

## License

MIT License. See [LICENSE](LICENSE).
