# NoHarm Workflow

NoHarm is designed as a low-friction triage layer for transcript isoform sets.
The first public workflow fixes the encoder so users do not need to tune a
model before seeing whether their data have an extreme isoform-divergence tail.

Version 0.1 is intentionally simple: it compares each isoform against a shared
uniform codon baseline, then compares the resulting structural residues within
each gene. It is not yet a gene-to-CDS or RNA-to-protein translation-alignment
model.

## Processing Diagram

```mermaid
flowchart TD
    A["Transcript isoforms"] --> B["Group isoforms<br/>by gene"]
    B --> C["Encode each isoform<br/>as 3-mer/codon windows"]
    C --> D["Compare with a shared baseline<br/>uniform 64-bin expectation"]
    D --> E["Read structural residue<br/>mean |codon-harm| per isoform"]
    E --> F["Compare residues within each gene<br/>max-min divergence"]
    F --> G["Rank candidate genes<br/>and contrast isoform pairs"]
```

## Minimal Command

```bash
noharm scan --fasta transcripts.fa.gz --out results/my_scan --workers 8
```

For a local, no-install run from the repository root:

```powershell
$env:PYTHONPATH="src"
python -m noharm scan --fasta data/demo_isoforms.fa --out results/demo
```

## Input Requirements

The default parser supports common transcript FASTA headers:

- GENCODE-style pipe-delimited headers where field 1 is transcript ID and field
  6 is gene name.
- `transcript|gene` two-field headers.
- headers containing `gene=...` or `gene_name=...`.

If a lab has a nonstandard FASTA format, the next practical extension is a
simple mapping table:

```text
transcript_id    gene_id
tx1              geneA
tx2              geneA
```

## Outputs

- `gene_divergence.tsv`: gene-level ranking.
- `isoform_scores.tsv`: per-transcript scores.
- `summary.json`: parameters, distribution, and top genes.
- `report.md`: compact report.

## Interpretation

The primary score asks:

> How different are the codon-landscape residues of isoforms from the same gene
> after projection onto a shared baseline?

High values nominate genes where isoform choice may produce a large codon-level
structural difference. Low values mean no large difference is detected in this
specific 3-mer coordinate; they do not prove functional equivalence. Translation
effects are a possible follow-up interpretation, not something directly modeled
by the v0.1 scanner.

## Release Caveats

- Current public scanner reproduces the primary `Delta |codon-harm|` ranking.
- It compares isoforms with a shared uniform baseline. It does not yet implement
  the deeper gene/transcript-to-CDS/protein alignment envisioned for later
  NoHarm work.
- Matched null and CDS/full controls are currently documented in
  `docs/EXPERIMENT_REPORT.md`; publication-grade biological claims require GC
  matching, GTF-based CDS extraction, and formal GO enrichment.

## Broader Theory

This repo is meant to be usable without accepting any broader theory. Readers
who want the conceptual background can browse the GBE project:
[https://jackeylgene.github.io/GBE](https://jackeylgene.github.io/GBE).
