"""Core NoHarm isoform divergence scanner."""

from __future__ import annotations

import csv
import json
import math
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from .encode import gc_fraction, norm, window_vectors
from .fasta import FastaRecord, read_fasta
from .frame import FrameEconomy


@dataclass(frozen=True)
class ScanParams:
    window_codons: int = 30
    min_codons: int = 5
    memory_cap: int = 16
    workers: int = 1
    batch_size: int = 500


@dataclass
class IsoformScore:
    gene: str
    transcript_id: str
    mean_ch: float
    tau: float
    length: int
    gc: float
    n_windows: int
    n_codons: int


@dataclass
class GeneScore:
    gene: str
    n_isoforms: int
    ch_range: float
    ch_std: float
    tau_range: float
    len_range: int
    gc_range: float
    mean_ch: float
    mean_tau: float
    min_ch_transcript: str
    max_ch_transcript: str
    min_ch: float
    max_ch: float


def score_record(record: FastaRecord, params: ScanParams) -> IsoformScore | None:
    vectors = window_vectors(record.sequence, params.window_codons, None, params.min_codons)
    if not vectors:
        return None
    mags = [norm(vec) for vec in vectors]
    frame = FrameEconomy(memory_cap=params.memory_cap)
    for vec in vectors:
        frame.process(vec)
    clean_len = len(record.sequence.replace("\n", "").replace("\r", ""))
    return IsoformScore(
        gene=record.gene,
        transcript_id=record.transcript_id,
        mean_ch=sum(mags) / len(mags),
        tau=frame.tau,
        length=clean_len,
        gc=gc_fraction(record.sequence),
        n_windows=len(vectors),
        n_codons=clean_len // 3,
    )


def _score_batch(batch: tuple[list[FastaRecord], ScanParams]) -> list[IsoformScore]:
    records, params = batch
    scores: list[IsoformScore] = []
    for record in records:
        score = score_record(record, params)
        if score is not None:
            scores.append(score)
    return scores


def _batches(items: list[FastaRecord], size: int) -> Iterable[list[FastaRecord]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def scan_fasta(fasta_path: str | Path, params: ScanParams) -> tuple[list[GeneScore], list[IsoformScore], dict]:
    t0 = time.time()
    records = list(read_fasta(fasta_path))
    total_records = len(records)

    if params.workers > 1:
        payloads = [(batch, params) for batch in _batches(records, params.batch_size)]
        isoforms: list[IsoformScore] = []
        with ProcessPoolExecutor(max_workers=params.workers) as executor:
            for result in executor.map(_score_batch, payloads):
                isoforms.extend(result)
    else:
        isoforms = []
        for record in records:
            score = score_record(record, params)
            if score is not None:
                isoforms.append(score)

    by_gene: dict[str, list[IsoformScore]] = defaultdict(list)
    for isoform in isoforms:
        by_gene[isoform.gene].append(isoform)

    genes: list[GeneScore] = []
    for gene, values in by_gene.items():
        if len(values) < 2:
            continue
        chs = [v.mean_ch for v in values]
        taus = [v.tau for v in values]
        lengths = [v.length for v in values]
        gcs = [v.gc for v in values]
        mean_ch = sum(chs) / len(chs)
        mean_tau = sum(taus) / len(taus)
        min_iso = min(values, key=lambda item: item.mean_ch)
        max_iso = max(values, key=lambda item: item.mean_ch)
        genes.append(
            GeneScore(
                gene=gene,
                n_isoforms=len(values),
                ch_range=max(chs) - min(chs),
                ch_std=math.sqrt(sum((v - mean_ch) ** 2 for v in chs) / len(chs)),
                tau_range=max(taus) - min(taus),
                len_range=max(lengths) - min(lengths),
                gc_range=max(gcs) - min(gcs),
                mean_ch=mean_ch,
                mean_tau=mean_tau,
                min_ch_transcript=min_iso.transcript_id,
                max_ch_transcript=max_iso.transcript_id,
                min_ch=min_iso.mean_ch,
                max_ch=max_iso.mean_ch,
            )
        )
    genes.sort(key=lambda item: item.ch_range, reverse=True)

    summary = build_summary(genes, isoforms, total_records, time.time() - t0, params)
    return genes, isoforms, summary


def percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return 0.0
    idx = min(len(sorted_values) - 1, int(q * len(sorted_values)))
    return sorted_values[idx]


def build_summary(
    genes: list[GeneScore],
    isoforms: list[IsoformScore],
    total_records: int,
    elapsed_seconds: float,
    params: ScanParams,
) -> dict:
    values = sorted(g.ch_range for g in genes)
    thresholds = {str(t): sum(v > t for v in values) for t in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]}
    return {
        "total_records": total_records,
        "scored_isoforms": len(isoforms),
        "multi_isoform_genes": len(genes),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "params": asdict(params),
        "ch_range_distribution": {
            "min": min(values) if values else 0.0,
            "p50": percentile(values, 0.50),
            "p90": percentile(values, 0.90),
            "p95": percentile(values, 0.95),
            "p99": percentile(values, 0.99),
            "max": max(values) if values else 0.0,
        },
        "threshold_counts": thresholds,
        "top_genes": [asdict(g) for g in genes[:30]],
    }


def write_outputs(genes: list[GeneScore], isoforms: list[IsoformScore], summary: dict, out_dir: str | Path) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "gene_divergence.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(genes[0]).keys()) if genes else list(GeneScore.__annotations__))
        writer.writeheader()
        for gene in genes:
            writer.writerow(asdict(gene))

    with (out / "isoform_scores.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(isoforms[0]).keys()) if isoforms else list(IsoformScore.__annotations__))
        writer.writeheader()
        for isoform in isoforms:
            writer.writerow(asdict(isoform))

    with (out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    write_report(out / "report.md", summary)


def write_report(path: Path, summary: dict) -> None:
    dist = summary["ch_range_distribution"]
    lines = [
        "# NoHarm Isoform Divergence Report",
        "",
        f"- Total FASTA records: {summary['total_records']}",
        f"- Scored isoforms: {summary['scored_isoforms']}",
        f"- Multi-isoform genes: {summary['multi_isoform_genes']}",
        f"- Elapsed seconds: {summary['elapsed_seconds']}",
        "",
        "## Distribution",
        "",
        f"- P50: {dist['p50']:.6f}",
        f"- P90: {dist['p90']:.6f}",
        f"- P95: {dist['p95']:.6f}",
        f"- P99: {dist['p99']:.6f}",
        f"- Max: {dist['max']:.6f}",
        "",
        "## Top Genes",
        "",
        "| Rank | Gene | Isoforms | Delta | Mean | Length Range | Top Isoform Pair |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for idx, gene in enumerate(summary["top_genes"][:30], 1):
        pair = f"{gene['min_ch_transcript']} -> {gene['max_ch_transcript']}"
        lines.append(
            f"| {idx} | {gene['gene']} | {gene['n_isoforms']} | "
            f"{gene['ch_range']:.6f} | {gene['mean_ch']:.6f} | "
            f"{gene['len_range']} | {pair} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "NoHarm ranks genes by the divergence among their transcript isoforms in a fixed 3-mer frame-economy coordinate.",
            "The score is an annotation-free prioritization signal, not a clinical or causal claim.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

