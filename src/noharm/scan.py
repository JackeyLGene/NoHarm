"""Core NoHarm isoform divergence scanner.

Ranking metrics:
  ch_range      - static L2 distance from uniform baseline
  tau_range     - frame-economy stability trace range
  merge_range   - frame-economy merge rate range
  novelty_range - frame-economy novelty stability range
  drift_range   - frame-economy centroid drift range
  churn_range   - frame-economy frame churn range
"""

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

# Metric keys exposed for ranking
FE_METRICS = ["ch_range", "tau_range", "merge_range", "novelty_range", "drift_range", "churn_range"]


@dataclass(frozen=True)
class ScanParams:
    window_codons: int = 30
    min_codons: int = 5
    memory_cap: int = 16
    merge_radius: float = 0.08
    workers: int = 1
    batch_size: int = 500
    rank_by: str = "ch_range"


@dataclass
class IsoformScore:
    gene: str
    transcript_id: str
    # static encoding (L2 baseline)
    mean_ch: float
    # frame-economy dynamic metrics
    tau: float
    merge_rate: float
    novelty_stability: float
    centroid_drift: float
    frame_churn: float
    final_occupancy: float
    # metadata
    length: int
    gc: float
    n_windows: int
    n_codons: int

    def get(self, key: str) -> float:
        return getattr(self, {"ch_range": "mean_ch", "tau_range": "tau",
            "merge_range": "merge_rate", "novelty_range": "novelty_stability",
            "drift_range": "centroid_drift", "churn_range": "frame_churn"}.get(key, key))


@dataclass
class GeneScore:
    gene: str
    n_isoforms: int
    ch_range: float
    ch_std: float
    tau_range: float
    merge_range: float
    novelty_range: float
    drift_range: float
    churn_range: float
    len_range: int
    gc_range: float
    mean_ch: float
    mean_tau: float
    min_ch_transcript: str
    max_ch_transcript: str
    min_ch: float
    max_ch: float
    contrast_metric: str
    min_contrast_transcript: str
    max_contrast_transcript: str
    min_contrast: float
    max_contrast: float


def score_record(record: FastaRecord, params: ScanParams) -> IsoformScore | None:
    vectors = window_vectors(record.sequence, params.window_codons, None, params.min_codons)
    if not vectors:
        return None
    mags = [norm(vec) for vec in vectors]
    frame = FrameEconomy(memory_cap=params.memory_cap, merge_radius=params.merge_radius)
    frame.warmup(params.memory_cap * 2)  # eliminate cold-start artifact
    for vec in vectors:
        frame.process(vec)
    clean_len = len(record.sequence.replace("\n", "").replace("\r", ""))
    return IsoformScore(
        gene=record.gene,
        transcript_id=record.transcript_id,
        mean_ch=sum(mags) / len(mags),
        tau=frame.tau,
        merge_rate=frame.merge_rate,
        novelty_stability=frame.novelty_stability,
        centroid_drift=frame.centroid_drift,
        frame_churn=frame.frame_churn,
        final_occupancy=frame.final_occupancy,
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


def _range_and_std(values: list[float]) -> tuple[float, float]:
    r = max(values) - min(values)
    mean = sum(values) / len(values)
    std = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    return r, std


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
        merges = [v.merge_rate for v in values]
        novelties = [v.novelty_stability for v in values]
        drifts = [v.centroid_drift for v in values]
        churns = [v.frame_churn for v in values]
        lengths = [v.length for v in values]
        gcs = [v.gc for v in values]

        ch_range, ch_std = _range_and_std(chs)
        tau_range, _ = _range_and_std(taus)
        merge_range, _ = _range_and_std(merges)
        novelty_range, _ = _range_and_std(novelties)
        drift_range, _ = _range_and_std(drifts)
        churn_range, _ = _range_and_std(churns)

        mean_ch = sum(chs) / len(chs)
        mean_tau = sum(taus) / len(taus)
        min_iso = min(values, key=lambda item: item.mean_ch)
        max_iso = max(values, key=lambda item: item.mean_ch)
        contrast_metric = params.rank_by
        min_contrast_iso = min(values, key=lambda item: item.get(contrast_metric))
        max_contrast_iso = max(values, key=lambda item: item.get(contrast_metric))

        genes.append(GeneScore(
            gene=gene, n_isoforms=len(values),
            ch_range=ch_range, ch_std=ch_std,
            tau_range=tau_range, merge_range=merge_range,
            novelty_range=novelty_range, drift_range=drift_range,
            churn_range=churn_range,
            len_range=max(lengths) - min(lengths),
            gc_range=max(gcs) - min(gcs),
            mean_ch=mean_ch, mean_tau=mean_tau,
            min_ch_transcript=min_iso.transcript_id,
            max_ch_transcript=max_iso.transcript_id,
            min_ch=min_iso.mean_ch, max_ch=max_iso.mean_ch,
            contrast_metric=contrast_metric,
            min_contrast_transcript=min_contrast_iso.transcript_id,
            max_contrast_transcript=max_contrast_iso.transcript_id,
            min_contrast=min_contrast_iso.get(contrast_metric),
            max_contrast=max_contrast_iso.get(contrast_metric),
        ))

    rank_key = params.rank_by
    genes.sort(key=lambda item: getattr(item, rank_key), reverse=True)

    summary = build_summary(genes, isoforms, total_records, time.time() - t0, params)
    return genes, isoforms, summary


def percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return 0.0
    idx = min(len(sorted_values) - 1, int(q * len(sorted_values)))
    return sorted_values[idx]


def _dist(values: list[float]) -> dict:
    s = sorted(values)
    return {
        "min": s[0] if s else 0.0,
        "p50": percentile(s, 0.50), "p90": percentile(s, 0.90),
        "p95": percentile(s, 0.95), "p99": percentile(s, 0.99),
        "max": s[-1] if s else 0.0,
    }


def build_summary(
    genes: list[GeneScore],
    isoforms: list[IsoformScore],
    total_records: int,
    elapsed_seconds: float,
    params: ScanParams,
) -> dict:
    ch_vals = [g.ch_range for g in genes]
    tau_vals = [g.tau_range for g in genes]
    rank_key = params.rank_by
    rank_vals = [getattr(g, rank_key) for g in genes]
    thresholds = {str(t): sum(v > t for v in ch_vals) for t in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]}
    dists = {k: _dist([getattr(g, k) for g in genes]) for k in FE_METRICS if hasattr(genes[0], k)} if genes else {}
    return {
        "total_records": total_records,
        "scored_isoforms": len(isoforms),
        "multi_isoform_genes": len(genes),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "params": asdict(params),
        "rank_by": rank_key,
        "distributions": dists,
        "threshold_counts": thresholds,
        "top_genes": [asdict(g) for g in genes[:30]],
    }


def write_outputs(genes: list[GeneScore], isoforms: list[IsoformScore], summary: dict, out_dir: str | Path) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "gene_divergence.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(genes[0]).keys()) if genes else [], delimiter="\t")
        writer.writeheader()
        for gene in genes:
            writer.writerow(asdict(gene))

    with (out / "isoform_scores.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(isoforms[0]).keys()) if isoforms else [], delimiter="\t")
        writer.writeheader()
        for isoform in isoforms:
            writer.writerow(asdict(isoform))

    with (out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    write_report(out / "report.md", summary)


def write_report(path: Path, summary: dict) -> None:
    rank_key = summary.get("rank_by", "ch_range")
    dist = summary["distributions"].get(rank_key, {})
    lines = [
        "# NoHarm Isoform Divergence Report",
        "",
        f"- Ranking metric: {rank_key}",
        f"- Total FASTA records: {summary['total_records']}",
        f"- Scored isoforms: {summary['scored_isoforms']}",
        f"- Multi-isoform genes: {summary['multi_isoform_genes']}",
        f"- Elapsed seconds: {summary['elapsed_seconds']}",
        "",
        "## Distribution",
        "",
        f"- P50: {dist.get('p50', 0):.6f}",
        f"- P90: {dist.get('p90', 0):.6f}",
        f"- P95: {dist.get('p95', 0):.6f}",
        f"- P99: {dist.get('p99', 0):.6f}",
        f"- Max: {dist.get('max', 0):.6f}",
        "",
        "## Top Genes",
        "",
        "| Rank | Gene | Isoforms | Delta | Mean ch | Length Range | Contrast Pair |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for idx, gene in enumerate(summary["top_genes"][:30], 1):
        pair = f"{gene['min_contrast_transcript']} -> {gene['max_contrast_transcript']}"
        lines.append(
            f"| {idx} | {gene['gene']} | {gene['n_isoforms']} | "
            f"{gene.get(rank_key, 0):.6f} | {gene.get('mean_ch', 0):.6f} | "
            f"{gene.get('len_range', 0)} | {pair} |"
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        "NoHarm ranks genes by the divergence among their transcript isoforms across",
        "a static codon-landscape coordinate and frame-economy response coordinates.",
        "Use ch_range for the static residue and merge_range for the primary",
        "frame-economy processing-response readout.",
        "",
        "The score is an annotation-free prioritization signal, not a clinical or causal claim.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
