"""Command line interface for NoHarm."""

from __future__ import annotations

import argparse
from pathlib import Path

from .scan import FE_METRICS, ScanParams, scan_fasta, write_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="noharm",
        description="Annotation-free structural screening of transcript isoform and region divergence.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    scan = sub.add_parser("scan", help="scan a transcript FASTA and rank isoform-divergent genes")
    scan.add_argument("--fasta", required=True, help="transcript FASTA path (.fa, .fasta, or .gz)")
    scan.add_argument("--out", default="results/noharm_scan", help="output directory")
    scan.add_argument("--window-codons", type=int, default=30, help="sliding window size in codons")
    scan.add_argument("--min-codons", type=int, default=5, help="minimum transcript length in codons")
    scan.add_argument("--memory-cap", type=int, default=16, help="frame trace memory cap")
    scan.add_argument("--merge-radius", type=float, default=0.08, help="frame merge radius")
    scan.add_argument("--workers", type=int, default=1, help="parallel workers")
    scan.add_argument("--batch-size", type=int, default=500, help="records per worker batch")
    scan.add_argument("--rank-by", default="ch_range", choices=FE_METRICS,
                      help="ranking metric: ch_range (static L2 baseline), tau_range, merge_range, "
                           "novelty_range, drift_range, churn_range")
    scan.add_argument("--top", type=int, default=10, help="number of top genes to print")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "scan":
        params = ScanParams(
            window_codons=args.window_codons,
            min_codons=args.min_codons,
            memory_cap=args.memory_cap,
            merge_radius=args.merge_radius,
            workers=max(1, args.workers),
            batch_size=max(1, args.batch_size),
            rank_by=args.rank_by,
        )
        genes, isoforms, summary = scan_fasta(args.fasta, params)
        write_outputs(genes, isoforms, summary, args.out)
        rank_key = args.rank_by

        label_map = {
            "ch_range": "static L2 (ablation baseline)",
            "tau_range": "frame-economy tau (stability)",
            "merge_range": "frame-economy merge rate",
            "novelty_range": "frame-economy novelty stability",
            "drift_range": "frame-economy centroid drift",
            "churn_range": "frame-economy frame churn",
        }

        print(f"NoHarm scan complete: {Path(args.out).resolve()}")
        print(f"  Ranked by: {rank_key} ({label_map.get(rank_key, '')})")
        print(f"  Records: {summary['total_records']}")
        print(f"  Scored isoforms: {summary['scored_isoforms']}")
        print(f"  Multi-isoform genes: {summary['multi_isoform_genes']}")
        dist = summary["distributions"].get(rank_key, {})
        print(
            f"  Distribution: P50={dist.get('p50', 0):.6f} "
            f"P95={dist.get('p95', 0):.6f} P99={dist.get('p99', 0):.6f} Max={dist.get('max', 0):.6f}"
        )
        print("\nTop genes:")
        for i, gene in enumerate(summary["top_genes"][: args.top], 1):
            print(
                f"  {i:>2}. {gene['gene']:<16} "
                f"delta={gene.get(rank_key, 0):.6f} "
                f"n={gene['n_isoforms']:<4} "
                f"pair={gene['min_contrast_transcript']}->{gene['max_contrast_transcript']}"
            )


if __name__ == "__main__":
    main()
