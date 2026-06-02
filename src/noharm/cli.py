"""Command line interface for NoHarm."""

from __future__ import annotations

import argparse
from pathlib import Path

from .scan import ScanParams, scan_fasta, write_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="noharm",
        description="Annotation-free prioritization of structurally divergent transcript isoforms.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    scan = sub.add_parser("scan", help="scan a transcript FASTA and rank isoform-divergent genes")
    scan.add_argument("--fasta", required=True, help="transcript FASTA path (.fa, .fasta, or .gz)")
    scan.add_argument("--out", default="results/noharm_scan", help="output directory")
    scan.add_argument("--window-codons", type=int, default=30, help="sliding window size in codons")
    scan.add_argument("--min-codons", type=int, default=5, help="minimum transcript length in codons")
    scan.add_argument("--memory-cap", type=int, default=16, help="frame trace memory cap")
    scan.add_argument("--workers", type=int, default=1, help="parallel workers")
    scan.add_argument("--batch-size", type=int, default=500, help="records per worker batch")
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
            workers=max(1, args.workers),
            batch_size=max(1, args.batch_size),
        )
        genes, isoforms, summary = scan_fasta(args.fasta, params)
        write_outputs(genes, isoforms, summary, args.out)

        print(f"NoHarm scan complete: {Path(args.out).resolve()}")
        print(f"  Records: {summary['total_records']}")
        print(f"  Scored isoforms: {summary['scored_isoforms']}")
        print(f"  Multi-isoform genes: {summary['multi_isoform_genes']}")
        dist = summary["ch_range_distribution"]
        print(
            f"  Distribution: P50={dist['p50']:.6f} "
            f"P95={dist['p95']:.6f} P99={dist['p99']:.6f} Max={dist['max']:.6f}"
        )
        print("\nTop genes:")
        for i, gene in enumerate(summary["top_genes"][: args.top], 1):
            print(
                f"  {i:>2}. {gene['gene']:<16} "
                f"delta={gene['ch_range']:.6f} "
                f"n={gene['n_isoforms']:<4} "
                f"pair={gene['min_ch_transcript']}->{gene['max_ch_transcript']}"
            )


if __name__ == "__main__":
    main()

