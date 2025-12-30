# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:48:27 2023

@author: Administrator
"""

from collections import OrderedDict
from pathlib import Path
import os

_DEFAULT_BED_DIR = Path(__file__).resolve().parents[1] / "data" / "genome"
file_dir = Path(os.environ.get("GENOME_BED_DIR", str(_DEFAULT_BED_DIR)))
color = [
    "lightgreen",
    "lightgreen",
    "pink",
    "lightblue",
    "pink",
    "lightblue",
    "pink",
    "lightblue",
    "pink",
    "lightblue",
    "pink",
    "lightblue",
    "pink",
    "lightblue",
]


def _load_gene_sites(bed_name: str):
    bed_path = file_dir / bed_name
    if not bed_path.exists():
        raise FileNotFoundError(
            f"Missing BED file: {bed_path}\n"
            "You can generate the BEDs that this repo expects with:\n"
            f"  python3 {Path(__file__).resolve().parents[0] / 'make_genome_beds.py'}\n"
        )

    out = []
    with bed_path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            gene = cols[3]
            start = int(cols[1])
            end = int(cols[2])
            out.append((gene, start, end))
    return out


gene_sites_sars = _load_gene_sites("NC_045512.2.bed")
gene_sites_pedv = _load_gene_sites("MK584552.1.bed")
gene_sites_pdcov_kv = _load_gene_sites("KY363867.1.bed")
gene_sites_pdcov_kt = _load_gene_sites("KT336560.1.bed")
gene_sites_ibv = _load_gene_sites("MK878536.1.bed")
gene_sites_sars1 = _load_gene_sites("NC_004718.3.bed")
gene_sites_mers = _load_gene_sites("NC_019843.3.bed")

dict_gene_sites_sars = OrderedDict()
for gene in gene_sites_sars:
    dict_gene_sites_sars[gene[0]] = (gene[1], gene[2])

dict_gene_sites_pedv = OrderedDict()
for gene in gene_sites_pedv:
    dict_gene_sites_pedv[gene[0]] = (gene[1], gene[2])

dict_gene_sites_pdcov = OrderedDict()
for gene in gene_sites_pdcov_kt:
    dict_gene_sites_pdcov[gene[0]] = (gene[1], gene[2])

dict_gene_sites_ibv = OrderedDict()
for gene in gene_sites_ibv:
    dict_gene_sites_ibv[gene[0]] = (gene[1], gene[2])

dict_gene_sites_sars1 = OrderedDict()
for gene in gene_sites_sars1:
    dict_gene_sites_sars1[gene[0]] = (gene[1], gene[2])

dict_gene_sites_mers = OrderedDict()
for gene in gene_sites_mers:
    dict_gene_sites_mers[gene[0]] = (gene[1], gene[2])

