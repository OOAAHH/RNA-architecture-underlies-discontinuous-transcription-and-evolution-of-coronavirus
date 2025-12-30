# -*- coding: utf-8 -*-
"""
Plot a long-range interaction figure (arc diagram + triangular heatmap) from a 3-column
interaction matrix file: `S  E  value` (coordinates in nt).

Example (SARS-CoV-2, 27,250–29,500 nt):
  MPLCONFIGDIR=.mplconfig python3 git-scripts/plot_long_range_interactions.py \
    --matrix "data/RNA-RNA interactions/SARS-CoV-2/Sars-cov2-genome-sample.5nt.kr.matrix" \
    --start 27250 --end 29500 --res 5 \
    --out figure/SARS2_27250_29500_long_range.png
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Interaction:
    s: int
    e: int
    value: float


def _iter_interactions(path: Path) -> Iterable[Interaction]:
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            s = int(float(parts[0]))
            e = int(float(parts[1]))
            v = float(parts[2])
            if s > e:
                s, e = e, s
            yield Interaction(s=s, e=e, value=v)


def _build_submatrix(
    interactions: Iterable[Interaction],
    start_nt: int,
    end_nt: int,
    res: int,
) -> Tuple[np.ndarray, List[Interaction], int, int]:
    start_bin = int(start_nt // res)
    end_bin = int(end_nt // res)
    bins = end_bin - start_bin + 1
    mat = np.zeros((bins, bins), dtype=np.float64)
    kept: List[Interaction] = []

    for it in interactions:
        if it.s < start_nt or it.e > end_nt:
            continue
        s_bin = int(round(it.s / res)) - start_bin
        e_bin = int(round(it.e / res)) - start_bin
        if s_bin < 0 or e_bin < 0 or s_bin >= bins or e_bin >= bins:
            continue
        mat[s_bin, e_bin] += it.value
        if s_bin != e_bin:
            mat[e_bin, s_bin] += it.value
        kept.append(it)

    return mat, kept, start_bin, end_bin


def _triangle_view(mat: np.ndarray) -> np.ndarray:
    l = mat.shape[0]
    h = l // 2
    tri = np.full((h, l), np.nan, dtype=np.float64)
    for d in range(h):
        row = h - d - 1
        for mid in range(d, l - d):
            i = mid - d
            j = mid + d
            tri[row, mid] = mat[i, j]
    return tri


def _plot_arc(ax, interactions: List[Interaction], min_distance: int, top_n: int) -> None:
    cand = [it for it in interactions if (it.e - it.s) >= min_distance and it.value > 0]
    cand.sort(key=lambda x: x.value, reverse=True)
    cand = cand[:top_n]

    for it in cand:
        s = it.s
        e = it.e
        radius = (e - s) / 2.0
        if radius <= 0:
            continue
        a = (s + e) / 2.0
        x = np.linspace(a - radius, a + radius, 256)
        y = np.sqrt(np.maximum(0.0, radius**2 - (x - a) ** 2))
        ax.plot(x, y, lw=0.25, color="grey", alpha=0.35)

    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)


def main(argv: Iterable[str] | None = None) -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_matrix = repo_root / "data" / "RNA-RNA interactions" / "SARS-CoV-2" / "Sars-cov2-genome-sample.5nt.kr.matrix"

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", type=Path, default=default_matrix, help="3-column interaction file (S E value).")
    p.add_argument("--start", type=int, default=27250, help="Region start (nt).")
    p.add_argument("--end", type=int, default=29500, help="Region end (nt).")
    p.add_argument("--res", type=int, default=5, help="Bin size (nt).")
    p.add_argument("--min-distance", type=int, default=2500, help="Minimum distance (nt) to draw arcs.")
    p.add_argument("--top-arcs", type=int, default=250, help="Draw at most N arcs.")
    p.add_argument("--cmap", type=str, default="viridis", help="Matplotlib colormap for heatmap.")
    p.add_argument("--vmax-percentile", type=float, default=99.0, help="Vmax percentile for log2(matrix+1).")
    p.add_argument(
        "--out",
        type=Path,
        default=repo_root / "figure" / "SARS2_27250_29500_long_range.png",
        help="Output image path.",
    )
    args = p.parse_args(list(argv) if argv is not None else None)

    if args.start >= args.end:
        raise SystemExit("--start must be < --end")
    if args.res <= 0:
        raise SystemExit("--res must be positive")
    if not args.matrix.exists():
        raise SystemExit(f"Matrix file not found: {args.matrix}")

    mat, kept, _start_bin, _end_bin = _build_submatrix(
        _iter_interactions(args.matrix),
        start_nt=args.start,
        end_nt=args.end,
        res=args.res,
    )

    tri = _triangle_view(np.log2(mat + 1.0))
    h, l = tri.shape

    nonzero = tri[np.isfinite(tri) & (tri > 0)]
    vmax = float(np.percentile(nonzero, args.vmax_percentile)) if nonzero.size else 1.0

    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.2], hspace=0.02)
    ax_arc = fig.add_subplot(gs[0])
    ax_tri = fig.add_subplot(gs[1], sharex=ax_arc)

    _plot_arc(ax_arc, kept, min_distance=args.min_distance, top_n=args.top_arcs)
    ax_arc.set_title(f"SARS-CoV-2: {args.start:,}-{args.end:,} nt", fontsize=14)

    extent = (args.start, args.end, 0, h * args.res)
    im = ax_tri.imshow(
        tri,
        origin="upper",
        aspect="auto",
        interpolation="nearest",
        cmap=plt.get_cmap(args.cmap),
        vmin=0,
        vmax=vmax,
        extent=extent,
    )
    ax_tri.set_ylim(h * args.res, 0)
    ax_tri.spines[["left", "right", "top"]].set_visible(False)
    ax_tri.set_yticks([])
    ax_tri.set_xlabel("Position (nt)")

    # Colorbar
    cax = fig.add_axes([0.92, 0.14, 0.02, 0.28])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log2(interaction + 1)")

    # Save
    out_path: Path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out_path}")
    return 0


if __name__ == "__main__":
    # Avoid matplotlib cache issues on locked-down machines (optional).
    os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))
    raise SystemExit(main())

