# RNA Architecture Underlies Discontinuous Transcription and Evolution of Coronavirus
RNA Architecture Underlies Discontinuous Transcription and Evolution of Coronavirus

Zi Wen#, Lei Chen#, Dehua Luo#, Ju Sun#, Liangrong Guo, Yingxiang Deng, Zhiyuan Huang, Yuxiang Wang, Ke Pan, Fan Wang*, Shaobo Xiao*, Li Li* and Dengguo Wei*

Correspondence: fwang@mail.ccnu.edu.cn; vet@mail.hzau.edu.cn; li.li@mail.hzau.edu.cn; dgwei@mail.hzau.edu.cn

Discontinuous transcription is employed by coronaviruses to generate subgenomic RNAs (sgRNAs) essential for gene expression and evolution adaptation. However, the mechanisms underlying the formation and functional roles of sgRNAs in specific coronavirus genomes are not well understood, particularly at the nucleotide level. We conducted an integrated analysis on 229 transcriptome and 12 RNA structurome samples across various coronavirus genera. RNA-RNA interactions between TRS-B and TRS-L flanking regions were identified in same genomic direction, which correlates with canonical junction formation. Non-canonical junctions frequently span either closely or distantly genomic distance, with short-range junctions generally mediated by stem-loops and overlapping with genomic deletion regions. Conserved long-range non-canonical sgRNAs were identified across different coronaviruses, and these sgRNAs harbor ORF10 or an evolving gene to suppress antiviral innate immune responses. This research highlights the importance of integrating transcriptome and RNA structurome profiles to elucidate the discontinuous transcription mechanism and evolution of coronaviruses.

Data can be  available on CODE (https://webofgroup.cn/dgwei/hngene).

## 我的新增
## Generate required `.bed` files

Some scripts in `git-scripts/` (for example `git-scripts/read_genome_sites.py`) expect genome annotation BED6 files under `data/genome/` (e.g. `NC_045512.2.bed`).

- Generate local BEDs (uses coordinate dictionaries already included in this repo): `python3 git-scripts/make_genome_beds.py`
- Generate the remaining BEDs by downloading GenBank records from NCBI (requires network access): `python3 git-scripts/make_genome_beds.py --ncbi`

## 新增了必要的用于复现的文件
我的目的是绘制SARS covid2长距离相互作用，图片已经复现，存在figure文件夹，必要的文件都存入仓库。
