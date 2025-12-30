# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:35:14 2023

@author: Administrator
"""
from Bio.Seq import Seq
from Bio import pairwise2 as pw2
from Bio.pairwise2 import format_alignment
from SARS2_recom_data import rna_seq_24_all
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
from ric_fun import classify_rna_data, get_max_line
from Bio.Align import substitution_matrices

sns.set_theme(style="whitegrid")
dir_pre = 'G:/' #'C:/Users/DELL/'  #
file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/Recobination/'
seq_sars = np.loadtxt(dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/NC_045512.2.fasta', skiprows=1, dtype = str, delimiter='/n')
seq_sars = np.array(list(''.join(seq_sars)))

recom_data = rna_seq_24_all.copy()
gap_lens = 100
rna_c, rna_nc, hist_data = classify_rna_data(recom_data, gap_lens, 10, 100, [100,2000])
distance = rna_nc['E'] - rna_nc['S']

l1 = 2500
l2 = 20000
virus_lens = 29903
l = 6
match = 1
mismatch = -1
open_tend = -1
extend = -0.5
name = "RNA"
matrix = substitution_matrices.load(name)

loop = rna_nc[(distance  < 2500)].copy()
rs = np.random.random_sample(rna_nc[(distance  < 2500)].shape[0]) * 27403
re = rs + np.random.random_sample(rna_nc[(distance  < 2500)].shape[0]) * 2500 + 100
rand = rna_nc[(distance  < 2500)].copy()
rand['S'] = rs
rand['E'] = re

list_score1 = [[],[],[]]
list_score2 = [[],[],[]]   
junc_value = []
files_out = open(file_dir + str(l1) + 'bp-loop-lf-rand-'+name+'-('+str(gap_lens)+'_'+str(match)+'_'+str(mismatch)+'_'+str(open_tend)+'_'+str(extend)+').txt','w')
for i in range(loop.shape[0]):
        start = int(loop['S'][i])
        end = int(loop['E'][i])
        start_r = int(rand['S'][i])
        end_r = int(rand['E'][i])
        value = loop['value'][i]
        values = int(loop['value'][i])
        if start < l or end > 29850:
            continue
        strs1 =  Seq(''.join(seq_sars[(start-l):(start + l)])) 
        strs2 =  Seq(''.join(seq_sars[(end-l):(end + l)]))
        strs_r1 =  Seq(''.join(seq_sars[(start_r-l):(start_r + l)])) 
        strs_r2 =  Seq(''.join(seq_sars[(end_r-l):(end_r + l)]))

        #match1 = pw2.align.localms(strs1,strs2,match, mismatch, open_tend,extend)[0]
        #match2 = pw2.align.localms(strs_r1,strs_r2,match, mismatch, open_tend,extend)[0]
        match1 = pw2.align.localds(strs1,strs2,matrix, open_tend,extend)[0]
        match2 = pw2.align.localds(strs_r1,strs_r2,matrix, open_tend,extend)[0]
        match_str1 = '\t'.join(list(map(str,match1)))
        match_str2 = '\t'.join(list(map(str,match2)))
        
        score1 = str(match1[2] / (match1[4] - match1[3]))
        score2 = str(match2[2] / (match2[4] - match2[3]))
        list_score1[0].append(float(score1))
        list_score1[1].append(match1[2])
        list_score2[0].append(float(score2))
        list_score2[1].append(match2[2])
        junc_value.append(value)
        
        paried1 = get_max_line(format_alignment(*match1).split('\n')[1])
        paried2 = get_max_line(format_alignment(*match2).split('\n')[1])
        list_score1[2].append(paried1)
        list_score2[2].append(paried2)
        
        files_out.writelines('>'+str(start)+'-'+str(end)+':'+str(value)+'\n')
        out_lines = strs1+'\t'+str(strs2)+'\t'+score1+'\t'+score2+'\t'+str(values)+'\n'
        files_out.writelines(out_lines)
        files_out.writelines(format_alignment(*match1))
        files_out.writelines(format_alignment(*match2)+'\n')

files_out.close()

labels = ['Discontinuous transcription','Random']
fig, ([ax1, ax2, ax3],[ax4, ax5, ax6]) = plt.subplots(nrows=2, ncols=3, figsize=(13, 8))
bplot1 = ax1.boxplot([list_score1[1],list_score2[1]],notch=True,  # notch shape
                     vert=True, patch_artist=True, labels=labels)
ax1.text(1.05, 8, 'P-value=%.2e'%ss.ranksums(list_score1[1], list_score2[1], alternative='greater')[1])
ax1.set_title('<2500bp')
ax1.set_ylabel('sequence consistency(bp)')

ax4.hist(list_score1[1], alpha = 0.5, label = 'junctions')
ax4.hist(list_score2[1], alpha = 0.5, label = labels[1])
ax4.legend()

print(np.mean(list_score1[1]),np.mean(list_score2[1]),np.median(list_score1[1]),np.median(list_score2[1]))

loop = rna_nc[(distance  < 20000) & (distance  > 2500)]
rs = np.random.random_sample(rna_nc[(distance  < 20000) & (distance  > 2500)].shape[0]) * 10000 
re = rs + np.random.random_sample(rna_nc[(distance  < 20000) & (distance  > 2500)].shape[0]) * 17500 + 2403
rand = rna_nc[(distance  < 20000) & (distance  > 2500)].copy()
rand['S'] = rs
rand['E'] = re

list_score1 = [[],[],[]]
list_score2 = [[],[],[]]   
junc_value = []
files_out = open(file_dir + str(l1)+'-'+str(l2)+'bp-loop-lf-rand-'+name+'-('+str(gap_lens)+'_'+str(match)+'_'+str(mismatch)+'_'+str(open_tend)+'_'+str(extend)+').txt','w')
l = 10
for i in range(loop.shape[0]):
        start = int(loop['S'][i])
        end = int(loop['E'][i])
        start_r = int(rand['S'][i])
        end_r = int(rand['E'][i])
        value = loop['value'][i]
        values = int(loop['value'][i])
        if start < l or end > 29850:
            continue
        if start_r < l :
            start_r += 10

        strs1 =  Seq(''.join(seq_sars[(start-l):(start + l)])) 
        strs2 =  Seq(''.join(seq_sars[(end-l):(end + l)]))
        strs_r1 =  Seq(''.join(seq_sars[(start_r-l):(start_r + l)])) 
        strs_r2 =  Seq(''.join(seq_sars[(end_r-l):(end_r + l)]))

        #match1 = pw2.align.localms(strs1,strs2,match, mismatch, open_tend,extend)[0]
        #match2 = pw2.align.localms(strs_r1,strs_r2,match, mismatch, open_tend,extend)[0]
        match1 = pw2.align.localds(strs1,strs2,matrix, open_tend,extend)[0]
        match2 = pw2.align.localds(strs_r1,strs_r2,matrix, open_tend,extend)[0]
        match_str1 = '\t'.join(list(map(str,match1)))
        match_str2 = '\t'.join(list(map(str,match2)))
        
        score1 = str(match1[2] / (match1[4] - match1[3]))
        score2 = str(match2[2] / (match2[4] - match2[3]))
        list_score1[0].append(float(score1))
        list_score1[1].append(match1[2])
        list_score2[0].append(float(score2))
        list_score2[1].append(match2[2])
        junc_value.append(value)
        
        paried1 = get_max_line(format_alignment(*match1).split('\n')[1])
        paried2 = get_max_line(format_alignment(*match2).split('\n')[1])
        list_score1[2].append(paried1)
        list_score2[2].append(paried2)
        
        files_out.writelines('>'+str(start)+'-'+str(end)+':'+str(value)+'\n')
        out_lines = strs1+'\t'+str(strs2)+'\t'+score1+'\t'+score2+'\t'+str(values)+'\n'
        files_out.writelines(out_lines)
        files_out.writelines(format_alignment(*match1))
        files_out.writelines(format_alignment(*match2)+'\n')

files_out.close()

bplot2 = ax2.boxplot([list_score1[1],list_score2[1]],notch=True,  # notch shape
                     vert=True, patch_artist=True, labels=labels)
ax2.text(1.05, 11, 'P-value=%.2e'%ss.ranksums(list_score1[1], list_score2[1], alternative='greater')[1])
ax2.set_title('2500-20000bp')

ax5.hist(list_score1[1], alpha = 0.5, label = 'junctions')
ax5.hist(list_score2[1], alpha = 0.5, label = labels[1])
ax5.legend()

print(np.mean(list_score1[1]),np.mean(list_score2[1]),np.median(list_score1[1]),np.median(list_score2[1]))

loop = rna_nc[distance  > 20000]
lens = rna_nc[distance  > 20000].shape[0]
rs = np.random.random_sample(lens * 5) * 9903 
re = rs + np.random.random_sample(lens * 5) * 9903 + 20000
mask = re < 29903
rand = rna_nc[distance  > 20000].copy()
rand['S'] = rs[mask][:lens]
rand['E'] = re[mask][:lens]

list_score1 = [[],[],[]]
list_score2 = [[],[],[]]   
junc_value = []
files_out = open(file_dir + str(l2)+'bp-loop-lf-rand-'+name+'-('+str(gap_lens)+'_'+str(match)+'_'+str(mismatch)+'_'+str(open_tend)+'_'+str(extend)+').txt','w')
l = 10
for i in range(loop.shape[0]):
        start = int(loop['S'][i])
        end = int(loop['E'][i])
        start_r = int(rand['S'][i])
        end_r = int(rand['E'][i])
        value = loop['value'][i]
        values = int(loop['value'][i])
        if start < l or end > 29850:
            continue
        strs1 =  Seq(''.join(seq_sars[(start-l):(start + l)])) 
        strs2 =  Seq(''.join(seq_sars[(end-l):(end + l)]))
        if start_r < l :
            start_r += 10
        strs_r1 =  Seq(''.join(seq_sars[(start_r-l):(start_r + l)])) 
        strs_r2 =  Seq(''.join(seq_sars[(end_r-l):(end_r + l)]))

        #match1 = pw2.align.localms(strs1,strs2,match, mismatch, open_tend,extend)[0]
        #match2 = pw2.align.localms(strs_r1,strs_r2,match, mismatch, open_tend,extend)[0]
        match1 = pw2.align.localds(strs1,strs2,matrix, open_tend,extend)[0]
        match2 = pw2.align.localds(strs_r1,strs_r2,matrix, open_tend,extend)[0]
        match_str1 = '\t'.join(list(map(str,match1)))
        match_str2 = '\t'.join(list(map(str,match2)))
        
        score1 = str(match1[2] / (match1[4] - match1[3]))
        score2 = str(match2[2] / (match2[4] - match2[3]))
        list_score1[0].append(float(score1))
        list_score1[1].append(match1[2])
        list_score2[0].append(float(score2))
        list_score2[1].append(match2[2])
        junc_value.append(value)
        
        paried1 = get_max_line(format_alignment(*match1).split('\n')[1])
        paried2 = get_max_line(format_alignment(*match2).split('\n')[1])
        list_score1[2].append(paried1)
        list_score2[2].append(paried2)
        
        files_out.writelines('>'+str(start)+'-'+str(end)+':'+str(value)+'\n')
        out_lines = strs1+'\t'+str(strs2)+'\t'+score1+'\t'+score2+'\t'+str(values)+'\n'
        files_out.writelines(out_lines)
        files_out.writelines(format_alignment(*match1))
        files_out.writelines(format_alignment(*match2)+'\n')

files_out.close()

bplot3 = ax3.boxplot([list_score1[1],list_score2[1]],notch=True,  # notch shape
                     vert=True, patch_artist=True, labels=labels)
ax3.text(1.05, 11, 'P-value=%.2e'%ss.ranksums(list_score1[1], list_score2[1], alternative='greater')[1])
ax3.set_title('>20000bp')

ax6.hist(list_score1[1], alpha = 0.5, label = 'junctions')
ax6.hist(list_score2[1], alpha = 0.5, label = labels[1])
ax6.legend()

print(np.mean(list_score1[1]),np.mean(list_score2[1]),np.median(list_score1[1]),np.median(list_score2[1]))

colors = ['pink', 'lightblue', 'lightgreen', 'lightyellow', 'grey', 'pink', 'lightblue', 'lightgreen', 'lightyellow']
for bplot in (bplot1, bplot2, bplot3):
    for patch, color in zip(bplot['boxes'], colors[:2]):
        patch.set_facecolor(color)


