# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:46:35 2023

@author: Administrator
"""

import numpy as np
from collections import Counter

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/'
rna_seq_24_all = np.loadtxt(file_dir + 'Recobination/SARS-cov2-all-rep-rna.1nt.none.matrix', usecols=[0,1,2]
                  , dtype=[('S','<i4'),('E','<i4'),('value','<f8')])


def read_nanopore(file_dir, file_name):
    for file in file_name:
        file_in = open(file_dir + file)
        out_list = []
        for line in file_in:
            list_line = line.split('\t')
            lens = len(list_line)
            n = int(lens / 2)
            if n == 2 :
                out_list.append((list_line[1], list_line[2]))
        
            elif n >= 3 :
                for i in np.arange(n-1):
                    out_list.append((list_line[i*2+1], list_line[i*2+2]))
        file_in.close()
    out_counts = Counter(out_list)
    return np.asarray([(key[0],key[1], out_counts[key]) for key in  out_counts], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])
    
def hist_rna_data(input_rna, dis, freq, binsize, discard):
    in_data = input_rna.copy()
    mask = (in_data['S'] <= discard[0]) & (in_data['value']>discard[1])
    data = in_data[~mask].copy()
    distance = data['E'] - data['S']
    rna_long = data[distance >= dis]
    rna_long_hf = rna_long[rna_long['value'] > freq]
    hist_data = []
    s = 0
    e = binsize
    lens = data['E'].max()
    while True:
        distance_hf = rna_long_hf['E'] - rna_long_hf['S']
        mask = (distance_hf >= s) & (distance_hf <= e)
        hist_data.append(rna_long_hf['value'][mask].sum())
        s = e
        e += binsize
        if e >= lens:
            break
    return hist_data
    
    
    
nanopore_24_1 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-1.bed'])
nanopore_24_2 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-2.bed'])
nanopore_24_3 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-3.bed'])

nanopore_24_10_2_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-1-2.bed', 'Vero_mapped_24-2-2.bed', 'Vero_mapped_24-3-2.bed'])
nanopore_24_10_3_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-1-3.bed', 'Vero_mapped_24-2-3.bed', 'Vero_mapped_24-3-3.bed'])

nanopore_24_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_24-1.bed', 'Vero_mapped_24-2.bed', 'Vero_mapped_24-3.bed'])

nanopore_48_1 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-1.bed'])
nanopore_48_2 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-2.bed'])
nanopore_48_3 = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-3.bed'])
nanopore_48_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-1.bed', 'Vero_mapped_48-2.bed', 'Vero_mapped_48-3.bed'])

nanopore_48_10_2_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-1-2.bed', 'Vero_mapped_48-2-2.bed', 'Vero_mapped_48-3-2.bed'])
nanopore_48_10_3_all = read_nanopore(file_dir+'nanopore/cp/', ['Vero_mapped_48-1-3.bed', 'Vero_mapped_48-2-3.bed', 'Vero_mapped_48-3-3.bed'])

nanopore_kim_24 = read_nanopore(file_dir+'nanopore/kim/', ['VeroInf24h.covleader.bed'])


