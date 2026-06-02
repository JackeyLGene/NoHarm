"""A small frame-economy trace used by the standalone NoHarm scanner.

The first public scanner ranks genes by isoform-level |cross-harm| divergence.
This light trace keeps a deterministic temporal summary without importing the
full research Geruon implementation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence


def euclidean(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(sum((x - y) * (x - y) for x, y in zip(a, b)))


@dataclass
class FrameEconomy:
    """Tiny centroid memory that produces a tau-like stability trace."""

    memory_cap: int = 16
    merge_radius: float = 0.08
    tau: float = 0.60
    frames: list[tuple[tuple[float, ...], float]] = field(default_factory=list)

    def process(self, vec: Sequence[float]) -> None:
        current = tuple(float(v) for v in vec)
        if not self.frames:
            novelty = 1.0
            self.frames.append((current, 1.0))
        else:
            distances = [euclidean(current, frame) for frame, _ in self.frames]
            best_idx = min(range(len(distances)), key=distances.__getitem__)
            best_dist = distances[best_idx]
            novelty = min(1.0, best_dist)
            if best_dist <= self.merge_radius:
                frame, weight = self.frames[best_idx]
                new_weight = weight + 1.0
                merged = tuple((frame[i] * weight + current[i]) / new_weight for i in range(len(current)))
                self.frames[best_idx] = (merged, new_weight)
            else:
                self.frames.append((current, 1.0))
                if len(self.frames) > self.memory_cap:
                    self.frames.sort(key=lambda item: item[1], reverse=True)
                    self.frames = self.frames[: self.memory_cap]

        target_tau = 0.55 + 0.25 * novelty
        self.tau = 0.95 * self.tau + 0.05 * target_tau

