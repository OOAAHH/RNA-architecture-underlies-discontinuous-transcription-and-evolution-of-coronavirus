# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 10:27:44 2023

@author: Administrator
"""

import numpy as np
from ric_fun import get_matrix_npz, get_oe
import matplotlib.pyplot as plt
import seaborn as sns
from SARS2_sgRNA import nanopore_24_10_2_all, nanopore_48_10_3_all, nanopore_24_10_3_all, nanopore_48_10_2_all

def extract_matrix(file_dir, file_name, res):
    data = np.loadtxt(file_dir + file_name, usecols=[0,1,2], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])
    data['S'] = data['S'] / res
    data['E'] = data['E'] / res
    mat = get_matrix_npz(data, data)
    oe, dis = get_oe(mat)
    return mat, dis

def compress_matrix(data, res):
    data = data.copy()
    data['S'] = data['S'] / res
    data['E'] = data['E'] / res
    mat = get_matrix_npz(data, data)
    oe, dis = get_oe(mat)
    return mat, dis

def DownSampling(rmat,ratio = 2):
    sampling_ratio = ratio
    m = np.matrix(rmat)

    all_sum = m.sum(dtype='float')
    m = m.astype(np.float64)
    idx_prob = np.divide(m, all_sum,out=np.zeros_like(m), where=all_sum != 0)
    idx_prob = np.asarray(idx_prob.reshape(
        (idx_prob.shape[0]*idx_prob.shape[1],)))
    idx_prob = np.squeeze(idx_prob)

    sample_number_counts = int(all_sum/(2*sampling_ratio))
    # 0 1 2 ... 8
    id_range = np.arange(m.shape[0]*m.shape[1])
    id_x = np.random.choice(
        id_range, size=sample_number_counts, replace=True, p=idx_prob)

    sample_m = np.zeros_like(m)
    for i in np.arange(sample_number_counts):
        x = int(id_x[i]/m.shape[0])
        y = int(id_x[i] % m.shape[0])
        sample_m[x, y] += 1.0
    sample_m = np.transpose(sample_m) + sample_m
    # print("after sample :",sample_m.sum())

    return np.array(sample_m)

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/hic/Splash/' #'C:/Users/DELL/OneDrive/vRic-seq/SARS-CoV2/hyb_mat/'
res = 5
cell_24, dis_cell_24 = extract_matrix(file_dir, 'Cell_5to3.5nt.VCRT.matrix', res)
cell_24_n, dis_cell_24_n = extract_matrix(file_dir, 'N-Cell.5nt.vcrt.matrix', res)
# cell_24_true = cell_24  - cell_24_n * 3
# cell_24_true[cell_24_true < 0] = 0

cell_48, dis_cell_48 = extract_matrix(file_dir, 'Lysate_5to3.5nt.VCRT.matrix', res)
cell_48_n, dis_cell_48_n = extract_matrix(file_dir, 'N-Late-cell.5nt.vcrt.matrix', res)
# cell_48_true = cell_48 - cell_48_n * 3
# cell_48_true[cell_48_true < 0] = 0

ratio = cell_48_n.sum() / cell_24_n.sum()
# cell_48_d = DownSampling(cell_48, cell_48_true.sum() / cell_24_true.sum())

nanopore_24_ma = compress_matrix(nanopore_24_10_2_all, res)[0]
nanopore_48_ma = compress_matrix(nanopore_48_10_2_all, res)[0] 
nanopore_24_ma_3 = compress_matrix(nanopore_24_10_3_all, res)[0]
nanopore_48_ma_3 = compress_matrix(nanopore_48_10_3_all, res)[0] 

plt.figure(figsize=(10,6))
vmax = 5
plt.subplot(131)
s1 = int(250 / 5)
s2 = int(3000 / 5)
e1 = int(29000 / 5)
e2 = int(29750 / 5)
n = 20
interp = 'kaiser'
ma_24 = cell_24.copy() 
ma_48 = cell_48.copy()
lens = cell_24.shape[0] 
plt.imshow(ma_24[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('24 hpi ' + str(ma_24[s1:s2,e1:e2].sum()) )
plt.subplot(132)
plt.imshow(ma_48[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('48 hpi '+str(ma_48[s1:s2,e1:e2].sum()) )
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)
plt.subplot(133)
plt.imshow(ma_48[s1:s2,e1:e2] - ma_24[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('48 hpi - 24 hpi ' + str(ma_48[s1:s2,e1:e2].sum() - ma_24[s1:s2,e1:e2].sum()) )
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)

sns.set_theme(style="ticks")
plt.figure(figsize=(10,5))
vmax = 5
k = 100
s = int(27300 / 5)
e = int(29550 / 5)
plt.subplot(131)
plt.imshow( (np.tril(ma_24[s:e,s:e], k = -k) + np.triu(ma_24[s:e,s:e], k = k) ), cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('24 hpi '+ str(np.triu(ma_24[s:e,s:e], k = k).sum()))
plt.subplot(132)
plt.imshow( (np.tril(ma_48[s:e,s:e], k = -k) + np.triu(ma_48[s:e,s:e], k = k) ), cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('48 hpi '+ str(np.triu(ma_48[s:e,s:e], k = k).sum()))
plt.subplot(133)
plt.imshow( ((np.tril(ma_48[s:e,s:e], k = -k) + np.triu(ma_48[s:e,s:e], k = k) + 1) - (np.tril(ma_24[s:e,s:e], k = -k) + np.triu(ma_24[s:e,s:e], k = k) + 1) ), cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('48 hpi - 24 hpi '+ str(np.triu(ma_48[s:e,s:e], k = k).sum() - np.triu(ma_24[s:e,s:e], k = k).sum()))


###sgRNA
plt.figure(figsize=(10,6))
vmax = 5
plt.subplot(131)
s1 = int(250 / 5)
s2 = int(3000 / 5)
e1 = int(29000 / 5)
e2 = int(29750 / 5)
n = 20
interp = 'kaiser'
ma_24 = nanopore_24_ma.copy()
ma_48 = nanopore_48_ma.copy()
lens = cell_24.shape[0] 
plt.imshow(ma_24[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('24 hpi ORF10-2')
plt.subplot(132)
plt.imshow(ma_48[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('48 hpi ORF10-2')
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)
plt.subplot(133)
plt.imshow(cell_24[s1:s2,e1:e2], interpolation = interp, cmap = plt.get_cmap('bwr'), vmax=vmax, vmin = -vmax)
plt.colorbar(location = 'bottom', fraction =0.01, pad = 0.075)
plt.title('24 hpi splash' )
plt.yticks(np.arange(0,s2-s1+1,res*n),np.arange(s1,s2+1,res*n)*res)
plt.xticks(np.arange(0,e2-e1+1,res*n),np.arange(e1,e2+1,res*n)*res)

sns.set_theme(style="ticks")
ma_24 = nanopore_24_ma_3.copy()
ma_48 = nanopore_48_ma_3.copy()
plt.figure(figsize=(10,5))
vmax = 4
k = 100
s = int(27300 / 5)
e = int(29550 / 5)
plt.subplot(131)
plt.imshow( np.log2(ma_24[s:e,s:e] + 1), cmap = plt.get_cmap('bwr'), vmax=vmax/2, vmin = -vmax/2)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('24 hpi '+ str(np.triu(ma_24[s:e,s:e], k = k).sum()))
plt.subplot(132)
plt.imshow( np.log2(ma_48[s:e,s:e] + 1), cmap = plt.get_cmap('bwr'), vmax=vmax/2, vmin = -vmax/2)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('48 hpi '+ str(np.triu(ma_48[s:e,s:e], k = k).sum()))
plt.subplot(133)
plt.imshow( np.log2(cell_24[s:e,s:e] + 1), cmap = plt.get_cmap('bwr'), vmax=vmax*2, vmin = -vmax*2)
plt.yticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.xticks(np.arange(0,e-s+1,res*n),np.arange(s,e+1,res*n)*res)
plt.colorbar(location = 'bottom', shrink = 0.25, pad = 0.075)
plt.title('24 hpi splash')


