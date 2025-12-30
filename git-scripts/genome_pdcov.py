# -*- coding: utf-8 -*-
"""
Minimal gene coordinate map for PDCoV (porcine deltacoronavirus).

This dictionary is used to generate `KT336560.1.bed` for scripts such as
`git-scripts/read_genome_sites.py`.
"""

from collections import OrderedDict

# Coordinates are kept exactly as in the original analysis scripts in this repo.
dict_gene = OrderedDict(
    [
        ("leader", (0, 100)),
        ("UTR_5", (101, 539)),
        ("ORF1a", (540, 11414)),
        ("ORF1b", (11414, 19342)),
        ("S", (19324, 22803)),
        ("E", (22797, 23048)),
        ("M", (23041, 23694)),
        ("NS6", (23694, 23978)),
        ("N", (23999, 25027)),
        ("NS7", (24093, 24695)),
        ("New1", (24194, 24659)),
        ("New2", (25040, 25240)),
        ("UTR_3", (25241, 25420)),
    ]
)

