"""FASTA parsing utilities."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


GENE_PATTERNS = [
    re.compile(r"(?:gene_name|gene|gene_symbol)=([^;\s|]+)"),
    re.compile(r"(?:GN)=([^;\s|]+)"),
]


@dataclass(frozen=True)
class FastaRecord:
    transcript_id: str
    gene: str
    header: str
    sequence: str


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_header(header: str) -> tuple[str, str]:
    """Return (transcript_id, gene) from common transcript FASTA headers."""
    text = header.strip()
    first = text.split()[0] if text else "unknown"

    for pattern in GENE_PATTERNS:
        match = pattern.search(text)
        if match:
            return first.split("|")[0], match.group(1)

    parts = first.split("|")
    if len(parts) >= 6:
        transcript_id = parts[0]
        gene = parts[5] or parts[1] or transcript_id
        return transcript_id, gene
    if len(parts) == 2:
        return parts[0], parts[1]
    if len(parts) > 2:
        return parts[0], parts[-1] or parts[0]

    transcript_id = first
    gene = transcript_id.split(".")[0]
    return transcript_id, gene


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    path = Path(path)
    header: str | None = None
    chunks: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    tid, gene = parse_header(header)
                    yield FastaRecord(tid, gene, header, "".join(chunks))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            tid, gene = parse_header(header)
            yield FastaRecord(tid, gene, header, "".join(chunks))

