"""Cell-Averaging CFAR ship detector for Sentinel-1 GRD.

CFAR (Constant False Alarm Rate) is the standard baseline for SAR ship
detection: each pixel is compared against a local clutter estimate derived
from a ring of "training" cells surrounding (but not touching) it. A guard
band between the test cell and training cells keeps target energy out of
the clutter estimate.

We use a fast 2D-summed-area-table implementation so the detector runs in
O(H·W) regardless of window size — practical for full GRD scenes.

The threshold is expressed as a multiplicative factor over the local mean,
which under the lognormal-clutter assumption corresponds to a stable false
alarm rate. See e.g. Crisp 2004 ("The state-of-the-art in ship detection in
SAR imagery"), used by the xView3 reference implementations.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Detection:
    row: int
    col: int
    snr: float
    bbox: tuple[int, int, int, int]  # (row_min, col_min, row_max, col_max), inclusive
    peak: float


def cfar_ca(
    image: np.ndarray,
    guard: int = 4,
    train: int = 12,
    pfa_factor: float = 3.0,
    min_pixels: int = 3,
) -> list[Detection]:
    """Cell-averaging CFAR.

    Args:
        image:       2D float32 SAR amplitude (e.g. Sentinel-1 GRD VV linear).
        guard:       half-width of the guard band (pixels). 0 = no guard.
        train:       half-width of the training window outside the guard.
        pfa_factor:  multiplicative threshold over the local clutter mean.
        min_pixels:  drop blobs below this footprint to filter speckle.

    Returns one Detection per connected-component blob exceeding the local
    threshold.
    """
    if image.ndim != 2:
        raise ValueError(f"expected 2D array, got shape {image.shape}")
    img = image.astype(np.float64, copy=False)

    win = train + guard
    pad = win
    padded = np.pad(img, pad, mode="reflect")
    # SAT with leading zero row/col so every window sum is one D - B - C + A
    # over half-open slices, regardless of where the window lands.
    sat = np.zeros((padded.shape[0] + 1, padded.shape[1] + 1), dtype=np.float64)
    sat[1:, 1:] = padded.cumsum(0).cumsum(1)

    H, W = img.shape

    def box_sum(half: int) -> np.ndarray:
        side = 2 * half + 1
        r0 = pad - half
        c0 = pad - half
        A = sat[r0 : r0 + H, c0 : c0 + W]
        B = sat[r0 + side : r0 + side + H, c0 : c0 + W]
        C = sat[r0 : r0 + H, c0 + side : c0 + side + W]
        D = sat[r0 + side : r0 + side + H, c0 + side : c0 + side + W]
        return D - B - C + A

    outer = box_sum(win)
    inner = box_sum(guard)
    outer_count = (2 * win + 1) ** 2
    inner_count = (2 * guard + 1) ** 2
    train_sum = outer - inner
    train_count = outer_count - inner_count
    clutter_mean = train_sum / train_count
    # Avoid divide-by-zero on flat water:
    clutter_mean = np.maximum(clutter_mean, 1e-12)
    threshold = pfa_factor * clutter_mean
    mask = img > threshold

    # Connected-component labeling (4-connectivity, pure numpy/scipy-free).
    detections: list[Detection] = []
    labels = _label_4(mask)
    n_labels = labels.max()
    if n_labels == 0:
        return detections
    snr_map = img / clutter_mean
    for lab in range(1, n_labels + 1):
        ys, xs = np.where(labels == lab)
        if ys.size < min_pixels:
            continue
        r0, r1 = int(ys.min()), int(ys.max())
        c0, c1 = int(xs.min()), int(xs.max())
        # Use the brightest pixel in the blob as the representative location.
        sub = img[r0 : r1 + 1, c0 : c1 + 1]
        peak_idx = np.unravel_index(np.argmax(sub), sub.shape)
        pr, pc = int(peak_idx[0]) + r0, int(peak_idx[1]) + c0
        detections.append(
            Detection(
                row=pr,
                col=pc,
                snr=float(snr_map[pr, pc]),
                bbox=(r0, c0, r1, c1),
                peak=float(img[pr, pc]),
            )
        )
    log.info("CFAR found %d detections (image %dx%d)", len(detections), H, W)
    return detections


def _label_4(mask: np.ndarray) -> np.ndarray:
    """4-connectivity connected components. Two-pass union-find. Pure numpy."""
    H, W = mask.shape
    labels = np.zeros((H, W), dtype=np.int32)
    parent: list[int] = [0]

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if ra < rb:
            parent[rb] = ra
        else:
            parent[ra] = rb

    next_label = 1
    for r in range(H):
        for c in range(W):
            if not mask[r, c]:
                continue
            up = labels[r - 1, c] if r > 0 else 0
            left = labels[r, c - 1] if c > 0 else 0
            if up == 0 and left == 0:
                labels[r, c] = next_label
                parent.append(next_label)
                next_label += 1
            elif up != 0 and left == 0:
                labels[r, c] = up
            elif up == 0 and left != 0:
                labels[r, c] = left
            else:
                lo = min(up, left)
                labels[r, c] = lo
                union(up, left)

    # Second pass: flatten labels via union-find roots, then renumber densely.
    remap: dict[int, int] = {}
    next_id = 0
    for i in range(1, next_label):
        root = find(i)
        if root not in remap:
            next_id += 1
            remap[root] = next_id
    flat = np.zeros_like(labels)
    for r in range(H):
        for c in range(W):
            v = labels[r, c]
            if v:
                flat[r, c] = remap[find(v)]
    return flat
