# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:39:39 2023

@author: Administrator
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from ric_fun import extract_matrix
#from plots import joint_heatmap_hist, joint_heatmap_hist_genome
from plots import joint_heatmap_genome
from read_genome_sites import dict_gene_sites_sars1, dict_gene_sites_ibv, dict_gene_sites_pdcov, dict_gene_sites_pedv, dict_gene_sites_sars, dict_gene_sites_mers
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

res = 50

# Default to this repository's data folder. Override with:
#   export JUNCTION_MATRIX_DIR="/path/to/data/juction sites"
_REPO_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_MATRIX_DIR = _REPO_ROOT / "data" / "juction sites"
file_dir = Path(os.environ.get("JUNCTION_MATRIX_DIR", str(_DEFAULT_MATRIX_DIR)))


def _dir(p: Path) -> str:
    # ric_fun.extract_matrix uses `file_dir + file_name`
    return str(p) + os.sep


def _pick_first_existing(subdir: str, candidates):
    base = file_dir / subdir
    for name in candidates:
        if (base / name).exists():
            return _dir(base), name
    raise FileNotFoundError(f"No matrix file found under: {base} (tried: {list(candidates)})")


# ma_pedv_rna_cell_1 = extract_matrix(_dir(file_dir / "PEDV"), "PEDV_rna-cell-1.1nt.none.matrix", res)[0]
ma_pedv_rna_cell = extract_matrix(_dir(file_dir / "PEDV"), "PEDV-rna-cell-2rep.1nt.none.matrix", res)[0]
ma_sars_rna = extract_matrix(_dir(file_dir / "SARS-CoV-2"), "SARS-cov2-all-rep-rna.1nt.none.matrix", res)[0]
ma_ibv_rna_cell = extract_matrix(_dir(file_dir / "IBV"), "IBV-rna.1nt.none.matrix", res)[0]
ma_sars1_rna = extract_matrix(_dir(file_dir / "SARS-CoV"), "Sars-cov-24hpi.1nt.none.matrix", res)[0]
ma_pdcov_rna = extract_matrix(_dir(file_dir / "PDCoV"), "pdcov-3rep.1nt.none.matrix", res)[0]
mers_dir, mers_name = _pick_first_existing(
    "MERS-CoV",
    [
        "MERS_SRR1942989.1nt.none.matrix",
        "MERS_SRR1195927.1nt.none.matrix",
    ],
)
ma_mers_rna = extract_matrix(mers_dir, mers_name, res)[0]

cmap = 'afmhot_r' #'PuBu'
# joint_heatmap_hist(ma_pedv_rna_cell, 'PEDV', 10, 5, cmap)
# joint_heatmap_hist(ma_sars_rna, 'SARS-CoV2', 10, 5, cmap)
# joint_heatmap_hist(ma_ibv_rna_cell, 'IBV', 10, 0.25, cmap)
# joint_heatmap_hist(ma_pdcov_rna, 'PDCoV', 10, 5, cmap)

# sns.set_theme(style="ticks")
# pp = PdfPages('G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/Fig-ORF10/geneme-junction-sites-new3.pdf')
# pp.savefig(joint_heatmap_hist_genome(ma_pdcov_rna, dict_gene_sites_pdcov, 'PDCoV', res, 10, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_ibv_rna_cell, dict_gene_sites_ibv, 'IBV', res, 2.5, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_pedv_rna_cell, dict_gene_sites_pedv, 'PEDV', res, 17.5, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_mers_rna, dict_gene_sites_mers, 'MERS-CoV', res, 15, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_sars_rna, dict_gene_sites_sars, 'SARS-CoV-2', res, 15, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_sars1_rna, dict_gene_sites_sars1, 'SARS-CoV', res, 8, cmap))
# pp.close()

pedv_utr3 = dict_gene_sites_pedv['UTR_3'][1] - dict_gene_sites_pedv['UTR_3'][0]
sars_utr3 = dict_gene_sites_sars['UTR_3'][1] - dict_gene_sites_sars['UTR_3'][0]
ibv_utr3 = dict_gene_sites_ibv['UTR_3'][1] - dict_gene_sites_ibv['UTR_3'][0]
pdcov_utr3 = dict_gene_sites_pdcov['UTR_3'][1] - dict_gene_sites_pdcov['UTR_3'][0]
sars1_utr3 = dict_gene_sites_sars1['UTR_3'][1] - dict_gene_sites_sars1['UTR_3'][0]
mers_utr3 = dict_gene_sites_mers['UTR_3'][1] - dict_gene_sites_mers['UTR_3'][0]

sns.set_theme(style="whitegrid")
fig = plt.figure(figsize=(8,6))
plt.bar(np.arange(1,7) - 0.155, [pedv_utr3, mers_utr3, sars_utr3, sars1_utr3, ibv_utr3, pdcov_utr3], .31, color = 'C0', alpha = 0.75, label = 'Length(bp)')
plt.xticks([1,2,3,4,5,6],['PEDV', 'MERS-CoV', 'SARS-CoV2', 'SARS-CoV', 'IBV', 'PDCoV'])
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.set_ylabel('Length(bp)')
ax2 = ax.twinx()
plt.bar(np.arange(1,7) + 0.155, [pedv_utr3/ dict_gene_sites_pedv['UTR_3'][1], mers_utr3/ dict_gene_sites_mers['UTR_3'][1], sars_utr3/ dict_gene_sites_sars['UTR_3'][1]
                    , sars1_utr3/ dict_gene_sites_sars['UTR_3'][1], ibv_utr3/  dict_gene_sites_ibv['UTR_3'][1], pdcov_utr3/  dict_gene_sites_pdcov['UTR_3'][1]], .31
        , color = 'C0', alpha = 0.5, label = 'Proportion')
ax2.spines['top'].set_visible(False)
ax2.set_ylabel('Proportion')
fig.legend(loc='upper center')

out_dir = Path(os.environ.get("FIGURE_DIR", str(_REPO_ROOT / "figure")))
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "junction_sites_utr3_length_proportion.png"
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved figure: {out_path}")


