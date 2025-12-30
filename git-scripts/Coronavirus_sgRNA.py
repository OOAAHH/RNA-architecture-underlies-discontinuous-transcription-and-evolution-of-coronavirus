# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 15:55:33 2023

@author: Administrator
"""

import numpy as np

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/'
rna_pedv_cell_all = np.loadtxt(file_dir + 'PEDV/matrix/PEDV-rna-cell-2rep.1nt.none.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_pdcov_cell_all = np.loadtxt(file_dir + 'PDCoV/matrix/WT_in_Virion.none.1bp.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_pdcov_cell_all_mut = np.loadtxt(file_dir + 'PDCoV/matrix/Mut_in_Virion.none.1bp.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

rna_seq_mers = np.loadtxt(file_dir + 'MERS-CoV/matrix/PRJNA279442-24h_MERS-rnaseq-PRJNA279442-calu3-24h.1nt.none.matrix'
                  , usecols=[0,1,2], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])

