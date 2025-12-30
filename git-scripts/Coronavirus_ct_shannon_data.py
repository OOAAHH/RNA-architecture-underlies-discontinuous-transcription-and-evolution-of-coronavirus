# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 20:31:41 2023

@author: Administrator
"""
import numpy as np
from Coronavirus_hic_data import ma_sars_splash_cell, ma_ric_pdcov_virion, oe_pedv_ric_virion_kr, oe_sars_ric_virion_kr, oe_ric_pdcov_virion
from ric_fun import get_shannon
res = 5
dir_pre = 'G:/'#'C:/Users/DELL/'

file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/CT/' 
strc_ric_sars2_virion = np.loadtxt( file_dir + 'virus.phase4_final.whole.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
shannon_ric_sars2_virion = np.loadtxt( file_dir + 'ric-shannon-sars-in-virion', usecols=[0], dtype='f8')


file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/PEDV/CT/' #'C:/Users/DELL/Desktop/vRic-seq/PEDV/'
strc_ric_pedv_cell = np.loadtxt( file_dir + 'PEDV-incell.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
strc_ric_pedv_virion = np.loadtxt( file_dir + 'PEDV-choice.phase4_final.whole.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)

# shannon_ric_pedv_cell = np.loadtxt( file_dir + 'ric-shannon-pedv-in-cell-KR-750', usecols=[0], dtype='f8')
# shannon_ric_pedv_virion = np.loadtxt( file_dir + 'ric-shannon-pedv-in-virion-KR-750', usecols=[0], dtype='f8')

file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/PDCOV/ct/' #'C:/Users/DELL/Desktop/vRic-seq/PEDV/'
strc_ric_pdcov_27_30 = np.loadtxt( file_dir + 'pdcov-27-30.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
strc_ric_pdcov_virion = np.loadtxt( file_dir + 'pdcov-ricseq-in-virion.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
strc_ric_pdcov_27_30_new = np.loadtxt( file_dir + 'PDCoV-in-virion-SciDat.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)


shannon_ric_pedv_virion = get_shannon(oe_pedv_ric_virion_kr, 150, 5)[1]
shannon_ric_sars2_virion = get_shannon(oe_sars_ric_virion_kr, 150, 5)[1]
shannon_ric_pdcov_virion = get_shannon(oe_ric_pdcov_virion, 150, 5)[1]


def get_connection(ma, strc, res):
    ma = ma.copy()
    strc_ric = strc.copy()
    value = ma[np.array((strc_ric['1']-1)/res, dtype = 'i4'),np.array((strc_ric['5']-1)/res, dtype = 'i4')]
    value[strc_ric['5'] == 0] = 0
    value = np.log(value / res**2 +1)
    return value

def save_connection(value, file_out):
    out = np.zeros(value.shape[0], dtype=[('Sites','<U6'),('Value','<f8')])
    out['Sites'] = np.arange(value.shape[0]) + 1
    out['Value'] = value
    return np.savetxt(file_out, out, delimiter='\t', fmt = '%s')


connection_splash = get_connection(ma_sars_splash_cell, strc_ric_sars2_virion, res)
connection_pdcov = get_connection(ma_ric_pdcov_virion, strc_ric_pdcov_virion, res)
