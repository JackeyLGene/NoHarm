"""Parity check against the original EE NoHarm genome scan.

This compares only the public ranking coordinate: per-isoform mean_ch and
per-gene ch_range. The standalone repo intentionally uses a lightweight tau
trace, so tau is not expected to match the original full Geruon.
"""

from __future__ import annotations

import gzip
import importlib.util
from collections import defaultdict
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from noharm.fasta import FastaRecord
from noharm.scan import ScanParams, score_record


EE_SCRIPT = Path(r"G:\GEME\EE\code\_noharm_genome_scan.py")
EE_CODE = EE_SCRIPT.parent
GENCODE_FASTA = Path(r"G:\GEME\EE\RNA\gencode\gencode.v49.pc_transcripts.fa.gz")


def load_original_module():
    sys.path.insert(0, str(EE_CODE))
    spec = importlib.util.spec_from_file_location("ee_noharm_genome_scan", EE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {EE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def first_records(limit: int = 200) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    with gzip.open(GENCODE_FASTA, "rt", encoding="utf-8") as handle:
        header = None
        chunks: list[str] = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    parts = header[1:].split("|")
                    if len(parts) >= 6:
                        records.append(FastaRecord(parts[0], parts[5], header[1:], "".join(chunks)))
                        if len(records) >= limit:
                            return records
                header = line
                chunks = []
            else:
                chunks.append(line)
        if header is not None and len(records) < limit:
            parts = header[1:].split("|")
            if len(parts) >= 6:
                records.append(FastaRecord(parts[0], parts[5], header[1:], "".join(chunks)))
    return records


def original_mean_ch(original, seq: str) -> float | None:
    n_codons = len(seq) // 3
    win = min(30, n_codons)
    if win < 5:
        return None
    stride = max(1, win // 5)
    mags = []
    for start in range(0, n_codons - win + 1, stride):
        begin = start * 3
        end = begin + win * 3
        actual = original.codon_freq_vec(seq[begin:end])
        hv = [actual[i] - original.UNIFORM_64[i] for i in range(64)]
        mags.append(sum(v * v for v in hv) ** 0.5)
    return sum(mags) / len(mags) if mags else None


def main() -> None:
    original = load_original_module()
    params = ScanParams()
    records = first_records()

    diffs = []
    by_gene_new = defaultdict(list)
    by_gene_old = defaultdict(list)
    compared = 0

    for record in records:
        old_ch = original_mean_ch(original, record.sequence)
        new = score_record(record, params)
        if old_ch is None or new is None:
            if old_ch is not None or new is not None:
                raise AssertionError(f"None mismatch for {record.transcript_id}: old={old_ch} new={new}")
            continue
        compared += 1
        new_ch = new.mean_ch
        diffs.append(abs(old_ch - new_ch))
        by_gene_old[record.gene].append(old_ch)
        by_gene_new[record.gene].append(new_ch)

    max_iso_diff = max(diffs) if diffs else 0.0
    if max_iso_diff > 1e-12:
        raise AssertionError(f"mean_ch mismatch: max diff={max_iso_diff}")

    gene_diffs = []
    for gene, old_values in by_gene_old.items():
        if len(old_values) < 2:
            continue
        new_values = by_gene_new[gene]
        old_range = max(old_values) - min(old_values)
        new_range = max(new_values) - min(new_values)
        gene_diffs.append(abs(old_range - new_range))

    max_gene_diff = max(gene_diffs) if gene_diffs else 0.0
    if max_gene_diff > 1e-12:
        raise AssertionError(f"ch_range mismatch: max diff={max_gene_diff}")

    print(f"parity passed: {compared} isoforms, max mean_ch diff={max_iso_diff:.3g}, max ch_range diff={max_gene_diff:.3g}")


if __name__ == "__main__":
    main()
