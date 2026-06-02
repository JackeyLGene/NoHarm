"""Sequence encoders for NoHarm.

The default encoder is intentionally fixed: non-overlapping 3-mer frequency
vectors in a sliding codon window, measured against a uniform 64-bin baseline.
This is the same low-prior signal used by the first full GENCODE scan.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence


BASE = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
UNIFORM_64 = tuple([1.0 / 64.0] * 64)


def normalize_sequence(seq: str) -> str:
    """Uppercase a DNA/RNA sequence while preserving frame positions.

    Ambiguous bases such as N are intentionally not removed. The original
    research scan skipped invalid 3-mers inside each window but did not collapse
    the sequence, so preserving positions is required for parity.
    """
    return seq.upper().replace("U", "T")


def clean_sequence(seq: str) -> str:
    """Return only canonical bases.

    This helper is kept for summary statistics. Encoding functions use
    ``normalize_sequence`` to preserve frame positions.
    """
    return "".join(ch for ch in normalize_sequence(seq) if ch in BASE)


def kmer3_index(codon: str) -> int:
    """Map a 3-mer over A/C/G/T to a stable 0..63 index."""
    return BASE[codon[0]] * 16 + BASE[codon[1]] * 4 + BASE[codon[2]]


def codon_freq_vec(seq: str) -> tuple[float, ...]:
    """Return a normalized 64-dimensional non-overlapping 3-mer vector."""
    vec = [0.0] * 64
    n = 0
    seq = normalize_sequence(seq)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3 and all(ch in BASE for ch in codon):
            vec[kmer3_index(codon)] += 1.0
            n += 1
    if n == 0:
        return tuple(vec)
    return tuple(v / n for v in vec)


def vector_delta(a: Sequence[float], b: Sequence[float]) -> tuple[float, ...]:
    return tuple(x - y for x, y in zip(a, b))


def norm(vec: Iterable[float]) -> float:
    return math.sqrt(sum(v * v for v in vec))


def window_vectors(
    seq: str,
    window_codons: int = 30,
    stride_codons: int | None = None,
    min_codons: int = 5,
) -> list[tuple[float, ...]]:
    """Encode a sequence as cross-harm window vectors.

    Each output vector is actual 3-mer frequency minus the uniform 64-bin
    baseline. For short transcripts, the window shrinks to the transcript
    length as long as at least ``min_codons`` are available.
    """
    seq = normalize_sequence(seq)
    n_codons = len(seq) // 3
    if n_codons < min_codons:
        return []
    win = min(window_codons, n_codons)
    stride = stride_codons if stride_codons is not None else max(1, win // 5)
    vectors: list[tuple[float, ...]] = []
    for start in range(0, n_codons - win + 1, stride):
        begin = start * 3
        end = begin + win * 3
        actual = codon_freq_vec(seq[begin:end])
        vectors.append(vector_delta(actual, UNIFORM_64))
    return vectors


def gc_fraction(seq: str) -> float:
    seq = clean_sequence(seq)
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)
