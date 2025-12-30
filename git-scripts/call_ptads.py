# -*- coding: utf-8 -*-
"""
Call pseudo-TADs (pTADs) from an RNA-RNA interaction map.

This implements a lightweight Directionality Index (DI)-style boundary calling
on a binned contact map:
  - For each bin i, compute upstream sum A and downstream sum B within a window.
  - Compute DI from A and B (Dixon et al. style).
  - Boundaries are called at robust DI sign changes; pTADs are segments between boundaries.

Supports:
  - `.mcool` (Cooler multi-resolution HDF5) with optional KR balancing (uses bins/KR).
  - 3-column `S E value` text files (like this repo's `*.matrix`).

Example:
  /path/to/python git-scripts/call_ptads.py \
    --mcool PRJNA655427/PRJNA655427_SARS2-vricseq-PRJNA655427-vero-NA-virion.mcool \
    --res 5 --window 1000 --norm KR --out figure/SARS2_virion_ptads_res5_win1000.bed
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class Ptad:
    start_nt: int
    end_nt: int


def _compute_di_from_pairs(
    pairs_bin1: np.ndarray,
    pairs_bin2: np.ndarray,
    values: np.ndarray,
    nbins: int,
    window_bins: int,
) -> np.ndarray:
    a = np.zeros(nbins, dtype=np.float64)
    b = np.zeros(nbins, dtype=np.float64)

    i = np.minimum(pairs_bin1, pairs_bin2)
    j = np.maximum(pairs_bin1, pairs_bin2)
    dist = j - i
    mask = (dist > 0) & (dist <= window_bins) & np.isfinite(values) & (values > 0)
    i = i[mask]
    j = j[mask]
    v = values[mask].astype(np.float64, copy=False)

    # For a symmetric matrix: for bin i, contacts to downstream (j) contribute to B[i];
    # for bin j, the same contacts are upstream and contribute to A[j].
    np.add.at(b, i, v)
    np.add.at(a, j, v)

    e = (a + b) / 2.0
    di = np.zeros(nbins, dtype=np.float64)
    nz = e > 0
    diff = b - a
    sign = np.sign(diff)
    # Dixon DI: sign * ((A-E)^2/E + (B-E)^2/E)
    di[nz] = sign[nz] * (((a[nz] - e[nz]) ** 2) / e[nz] + ((b[nz] - e[nz]) ** 2) / e[nz])
    di[~np.isfinite(di)] = 0.0
    return di


def _smooth(x: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return x
    w = int(window)
    kernel = np.ones(w, dtype=np.float64) / float(w)
    return np.convolve(x, kernel, mode="same")


def _call_boundaries(di: np.ndarray, min_abs_di: float, smooth_bins: int) -> np.ndarray:
    di2 = _smooth(di, smooth_bins)
    s = np.sign(di2)
    s[np.abs(di2) < min_abs_di] = 0

    boundaries = []
    prev = 0
    for idx in range(len(s)):
        cur = int(s[idx])
        if cur == 0:
            continue
        if prev != 0 and cur != prev:
            boundaries.append(idx)
        prev = cur
    return np.asarray(boundaries, dtype=np.int64)


def _ptads_from_boundaries(
    boundaries: np.ndarray,
    nbins: int,
    res: int,
    min_size_bins: int,
) -> List[Ptad]:
    raw_cuts = [0] + [int(x) for x in boundaries.tolist() if 0 < x < nbins] + [nbins]
    raw_cuts = sorted(set(raw_cuts))

    # Merge overly-fragmented boundaries: keep only cuts that produce segments >= min_size_bins.
    merged: List[Tuple[int, int]] = []
    seg_start = raw_cuts[0]
    for cut in raw_cuts[1:]:
        if cut - seg_start < min_size_bins:
            continue
        merged.append((seg_start, cut))
        seg_start = cut

    # Handle tail: merge into last segment if too short.
    if seg_start < nbins:
        if merged and (nbins - seg_start) < min_size_bins:
            merged[-1] = (merged[-1][0], nbins)
        elif (nbins - seg_start) >= min_size_bins:
            merged.append((seg_start, nbins))

    return [Ptad(start_nt=s * res, end_nt=e * res) for (s, e) in merged]


def _load_pairs_from_matrix(path: Path, res: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    # file format: S E value (S/E in nt; often multiples of res)
    s_list = []
    e_list = []
    v_list = []
    max_bin = 0
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
            bs = int(round(s / res))
            be = int(round(e / res))
            s_list.append(bs)
            e_list.append(be)
            v_list.append(v)
            max_bin = max(max_bin, bs, be)
    return (
        np.asarray(s_list, dtype=np.int64),
        np.asarray(e_list, dtype=np.int64),
        np.asarray(v_list, dtype=np.float64),
        max_bin + 1,
    )


def _load_pairs_from_mcool(path: Path, res: int, norm: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, str]:
    import h5py  # local import; present in the user's env

    with h5py.File(path, "r") as f:
        g = f["resolutions"][str(res)]
        chroms = g["chroms"]["name"][:]
        chrom = (
            chroms[0].decode("utf-8")
            if chroms.shape[0] > 0 and isinstance(chroms[0], (bytes, bytearray))
            else str(chroms[0]) if chroms.shape[0] > 0 else "HS1"
        )

        nbins = int(g["bins"]["start"].shape[0])
        b1 = g["pixels"]["bin1_id"][:].astype(np.int64, copy=False)
        b2 = g["pixels"]["bin2_id"][:].astype(np.int64, copy=False)
        c = g["pixels"]["count"][:].astype(np.float64, copy=False)

        if norm.upper() == "KR":
            w = g["bins"]["KR"][:].astype(np.float64, copy=False)
            v = c * w[b1] * w[b2]
        elif norm.lower() in {"none", "raw"}:
            v = c
        else:
            raise SystemExit(f"Unsupported --norm: {norm} (use KR or none)")

        v[~np.isfinite(v)] = 0.0
        return b1, b2, v, nbins, chrom


def main(argv: Optional[Iterable[str]] = None) -> int:
    repo_root = Path(__file__).resolve().parents[1]

    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--mcool", type=Path, help="Input .mcool file.")
    src.add_argument("--matrix", type=Path, help="Input 3-column S E value text file.")
    p.add_argument("--res", type=int, default=5, help="Resolution (nt per bin).")
    p.add_argument("--window", type=int, default=1000, help="DI window size in nt.")
    p.add_argument("--norm", type=str, default="KR", help="Normalization for .mcool: KR or none.")
    p.add_argument("--smooth", type=int, default=5, help="Smoothing window in bins for DI sign calling.")
    p.add_argument("--min-abs-di", type=float, default=0.0, help="Ignore DI with abs(di) below this.")
    p.add_argument("--min-size", type=int, default=500, help="Minimum pTAD size in nt.")
    p.add_argument("--out", type=Path, default=repo_root / "figure" / "ptads.bed", help="Output BED file.")
    args = p.parse_args(list(argv) if argv is not None else None)

    res = int(args.res)
    window_bins = max(1, int(round(args.window / res)))
    min_size_bins = max(1, int(round(args.min_size / res)))

    chrom = "HS1"
    if args.mcool:
        b1, b2, v, nbins, chrom = _load_pairs_from_mcool(args.mcool, res=res, norm=args.norm)
    else:
        b1, b2, v, nbins = _load_pairs_from_matrix(args.matrix, res=res)

    di = _compute_di_from_pairs(b1, b2, v, nbins=nbins, window_bins=window_bins)
    boundaries = _call_boundaries(di, min_abs_di=float(args.min_abs_di), smooth_bins=int(args.smooth))
    ptads = _ptads_from_boundaries(boundaries, nbins=nbins, res=res, min_size_bins=min_size_bins)

    out_path: Path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("chrom\tstart\tend\tname\n")
        for i, t in enumerate(ptads, start=1):
            f.write(f"{chrom}\t{t.start_nt}\t{t.end_nt}\tpTAD_{i}\n")

    print(f"Wrote {len(ptads)} pTADs to {out_path}")
    print(f"Boundaries (bins): {boundaries.tolist()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
