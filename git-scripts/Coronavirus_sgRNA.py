# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 15:55:33 2023

@author: Administrator
"""

import numpy as np

file_dir = '/Users/joseperezmartinez/docs/RNA-architecture-underlies-discontinuous-transcription-and-evolution-of-coronavirus-main/data'
rna_pedv_cell_all = np.loadtxt('/Users/joseperezmartinez/docs/RNA-architecture-underlies-discontinuous-transcription-and-evolution-of-coronavirus-main/data/juction sites/PEDV/PEDV-rna-cell-2rep.1nt.none.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_pdcov_cell_all = np.loadtxt('/Users/joseperezmartinez/docs/RNA-architecture-underlies-discontinuous-transcription-and-evolution-of-coronavirus-main/data/juction sites/PEDV/WT_in_Virion.none.1bp.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_pdcov_cell_all_mut = np.loadtxt('/Users/joseperezmartinez/docs/RNA-architecture-underlies-discontinuous-transcription-and-evolution-of-coronavirus-main/data/juction sites/PEDV/Mut_in_Virion.none.1bp.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_seq_mers = np.loadtxt('/Users/joseperezmartinez/docs/RNA-architecture-underlies-discontinuous-transcription-and-evolution-of-coronavirus-main/data/juction sites/PEDV/PRJNA279442-24h_MERS-rnaseq-PRJNA279442-calu3-24h.1nt.none.matrix'
                  , usecols=[0,1,2], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

