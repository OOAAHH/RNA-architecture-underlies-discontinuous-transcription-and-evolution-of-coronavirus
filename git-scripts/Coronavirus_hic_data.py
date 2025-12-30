# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 20:25:02 2023

@author: Administrator
"""
from ric_fun import extract_matrix
file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/Recobination/' #'C:/Users/61918/Desktop/'
file_dir_com = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/hic/Comrades/' #'C:/Users/61918/Desktop/'
file_dir_ric = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/hic/vRic-seq/' #'C:/Users/61918/Desktop/'
file_dir_splash = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/matrix/' #'C:/Users/61918/Desktop/'

res = 5
ma_sars_rna_cell, oe_sars_rna_cell, dis_sars_rna_cell = extract_matrix(file_dir, 'SARS-cov2-all-rep-rna.1nt.none.matrix', res)
ma_sars_com_cell, oe_sars_com_cell, dis_sars_com_cell = extract_matrix(file_dir_com, 'SARS-cov2-sample-minus.1nt.none.matrix', res)
ma_sars_ric_virion, oe_sars_ric_virion, dis_sars_ric_virion = extract_matrix(file_dir_ric, 'Sars-cov2-vRic-z.1nt.none.matrix', res)
ma_sars_ric_virion_kr, oe_sars_ric_virion_kr, dis_sars_ric_virion_kr = extract_matrix(file_dir_ric, 'SARS2_in_virion.5nt.KR.matrix', res)
ma_sars_splash_cell, oe_sars_splash_cell, dis_sars_splash_cell = extract_matrix(file_dir_splash, 'splash-SARS2-wanyue-WT.5nt.KR.matrix', res)
ma_sars_splash_382, oe_sars_splash_382, dis_sars_splash_382 = extract_matrix(file_dir_splash, 'splash-SARS2-wanyue-382.5nt.KR.matrix', res)

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/PDCoV/matrix/' #'C:/Users/61918/Desktop/'
res = 5
ma_ric_pdcov_cell, oe_ric_pdcov_cell, dis_ric_pdcov_cell = extract_matrix(file_dir, 'pdcov_in_Cell.5nt.KR.matrix', res)
ma_ric_pdcov_virion, oe_ric_pdcov_virion, dis_ric_pdcov_virion = extract_matrix(file_dir, 'PDCoV-27-30_in_Virion.5nt.kr.matrix', res)
ma_rna_pdcov_cell, oe_rna_pdcov_cell, dis_rna_pdcov_cell = extract_matrix(file_dir, 'WT_in_Virion.none.1bp.matrix', res)
ma_rna_pdcov_cell_mut, oe_rna_pdcov_cell_mut, dis_rna_pdcov_cell_mut = extract_matrix(file_dir, 'Mut_in_Virion.none.1bp.matrix', res)

file_dir_mers = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/MERS-CoV/matrix/' #'C:/Users/61918/Desktop/'
res = 5
ma_mers_com_cell_s, oe_mers_com_cell_s, dis_mers_com_cell_s = extract_matrix(file_dir_mers, 'MERS-ziv-Sample_in_Virion.5nt.KR.matrix', res)
ma_mers_com_cell_c, oe_mers_com_cell_c, dis_mers_com_cell_c = extract_matrix(file_dir_mers, 'MERS-ziv-Contrl_in_Virion.5nt.KR.matrix', res)
ma_mers_rna_cell, oe_mers_rna_cell, dis_mers_rna_cell = extract_matrix(file_dir_mers, 'PRJNA279442-24h_MERS-rnaseq-PRJNA279442-calu3-24h.1nt.none.matrix', res)
ma_mers_com_cell = ma_mers_com_cell_s - ma_mers_com_cell_c
ma_mers_com_cell[ma_mers_com_cell<0] = 0   
oe_mers_com_cell = oe_mers_com_cell_s - oe_mers_com_cell_c
oe_mers_com_cell[oe_mers_com_cell<0] = 0   

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/PEDV/matrix/' #'C:/Users/61918/Desktop/'
res = 5
ma_pedv_ric_cell_kr, oe_pedv_ric_cell_kr, dis_pedv_ric_cell_kr = extract_matrix(file_dir, 'PEDV_in_cell.5nt.KR.matrix', res)
ma_pedv_ric_cell_VCRT, oe_pedv_ric_cell_VCRT, dis_pedv_ric_cell_VCRT = extract_matrix(file_dir, 'PEDV_in_cell.5nt.VCRT.matrix', res)
ma_pedv_ric_cell, oe_pedv_ric_cell, dis_pedv_ric_cell = extract_matrix(file_dir, 'PEDV_in_cell-none-5nt.matrix', res)

ma_pedv_ric_virion_VCRT, oe_pedv_ric_virion_VCRT, dis_pedv_ric_virion_VCRT = extract_matrix(file_dir, 'PEDV_in_virion.5nt.VCRT.matrix', res)
ma_pedv_ric_virion_kr, oe_pedv_ric_virion_kr, dis_pedv_ric_virion_kr = extract_matrix(file_dir, 'PEDV_in_virion.5nt.KR.matrix', res)
ma_pedv_ric_virion, oe_pedv_ric_virion, dis_pedv_ric_virion = extract_matrix(file_dir, 'PEDV_in_virion.5nt.none.matrix', res)

ma_pedv_rna_cell, oe_pedv_rna_cell, dis_pedv_rna_cell = extract_matrix(file_dir, 'PEDV-rna-cell-2rep.1nt.none.matrix', res)
