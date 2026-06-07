"""Frame economy - finite centroid memory with temporal trace.

The FrameEconomy processes a sequence of window vectors through a
capacity-limited centroid memory. Each incoming vector either merges
into the nearest existing centroid (if within merge radius) or creates
a new frame (ejecting the weakest if at capacity).

Exposed metrics capture the memory's *dynamic response* to the stream,
not its static composition. These are the frame-economy readouts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence


def euclidean(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(sum((x - y) * (x - y) for x, y in zip(a, b)))


@dataclass
class FrameEconomy:
    """Finite centroid memory with configurable merge radius and capacity."""

    memory_cap: int = 16
    merge_radius: float = 0.08
    tau: float = 0.60

    # --- frame store ---
    frames: list[tuple[tuple[float, ...], float]] = field(default_factory=list)

    # --- dynamic trace ---
    total_windows: int = 0
    merge_count: int = 0
    novelty_sum: float = 0.0
    novelty_sq_sum: float = 0.0
    centroid_norms: list[float] = field(default_factory=list)
    frame_counts: list[int] = field(default_factory=list)

    def warmup(self, n_warm: int = 32) -> None:
        """Pre-feed uniform-zero vectors to fill memory and reach saturation.
        Eliminates cold-start artifact in per-isoform processing."""
        zero = tuple(0.0 for _ in range(64))
        for _ in range(n_warm):
            self._process_internal(zero)

    def process(self, vec: Sequence[float]) -> None:
        self._process_internal(tuple(float(v) for v in vec))

    def _process_internal(self, current: tuple) -> None:
        self.total_windows += 1

        if not self.frames:
            novelty = 1.0
            self.frames.append((current, 1.0))
        else:
            distances = [euclidean(current, f) for f, _ in self.frames]
            best_idx = min(range(len(distances)), key=distances.__getitem__)
            best_dist = distances[best_idx]
            novelty = min(1.0, best_dist / max(0.001, self.merge_radius))

            if best_dist <= self.merge_radius:
                self.merge_count += 1
                frame, weight = self.frames[best_idx]
                new_weight = weight + 1.0
                merged = tuple(
                    (frame[i] * weight + current[i]) / new_weight
                    for i in range(len(current))
                )
                self.frames[best_idx] = (merged, new_weight)
            else:
                self.frames.append((current, 1.0))
                if len(self.frames) > self.memory_cap:
                    self.frames.sort(key=lambda item: item[1], reverse=True)
                    self.frames = self.frames[: self.memory_cap]

        self.novelty_sum += novelty
        self.novelty_sq_sum += novelty * novelty
        self.centroid_norms.append(self._centroid_l2())
        self.frame_counts.append(len(self.frames))

        target_tau = 0.55 + 0.25 * novelty
        self.tau = 0.95 * self.tau + 0.05 * target_tau

    def _centroid_l2(self) -> float:
        if not self.frames:
            return 0.0
        total_weight = sum(w for _, w in self.frames)
        c = [0.0] * len(self.frames[0][0])
        for vec, w in self.frames:
            for i, v in enumerate(vec):
                c[i] += v * w
        for i in range(len(c)):
            c[i] /= max(total_weight, 1e-12)
        return math.sqrt(sum(x * x for x in c))

    # --- aggregate metrics ---

    @property
    def merge_rate(self) -> float:
        """Fraction of windows that merged into existing frames."""
        if self.total_windows == 0:
            return 0.0
        return self.merge_count / self.total_windows

    @property
    def novelty_stability(self) -> float:
        """Coefficient of variation of novelty. Lower = more stable processing."""
        if self.total_windows == 0:
            return 1.0
        mean = self.novelty_sum / self.total_windows
        variance = self.novelty_sq_sum / self.total_windows - mean * mean
        if variance < 0:
            variance = 0.0
        return math.sqrt(variance) / max(0.001, mean)

    @property
    def centroid_drift(self) -> float:
        """L2 distance between first and last centroid norm.
        Measures how far the memory center moved during processing."""
        if len(self.centroid_norms) < 2:
            return 0.0
        return abs(self.centroid_norms[-1] - self.centroid_norms[0])

    @property
    def frame_churn(self) -> float:
        """Std of frame count over time / mean frame count.
        Higher = more creation/pruning cycles."""
        if len(self.frame_counts) < 2:
            return 0.0
        mean = sum(self.frame_counts) / len(self.frame_counts)
        if mean == 0:
            return 0.0
        variance = sum((c - mean) ** 2 for c in self.frame_counts) / len(self.frame_counts)
        return math.sqrt(variance) / mean

    @property
    def final_occupancy(self) -> float:
        """Fraction of memory capacity currently occupied."""
        return len(self.frames) / max(1, self.memory_cap)
