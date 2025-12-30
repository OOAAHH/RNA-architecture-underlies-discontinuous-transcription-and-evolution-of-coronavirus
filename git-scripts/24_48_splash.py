# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:20:21 2024

@author: Administrator
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import OrderedDict
import scipy.stats as ss
from Corona_mRNA import run
sns.set_theme(style="whitegrid")

chrs = 'NC_045512'
file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/SARS-COV2/nanopore/'
dict_gene = OrderedDict([('leader', (0, 100)),
    ('UTR_5', (101, 265)),
    ('ORF1a', (266, 13468)),
    ('ORF1b', (13468, 21555)),
    ('S', (21563, 25384)),
    ('ORF3a', (25393, 26220)),
    ('E', (26245, 26492)),
    ('M', (26503, 27191)),
    ('ORF6', (27202, 27387)),
    ('ORF7a', (27394, 27759)),
    ('ORF7b', (27756, 27887)),
    ('ORF8', (27894, 28259)),
    ('N', (28274, 29533)),
    # ('ORF9b', (28284, 28577)),
    # ('ORF9c', (28734, 28955)),	
    ('ORF10', (29558, 29674)),
    ('UTR_3', (29675, 29903))])

tpm_c_rep1, reads_c_rep1, tpm_nc_rep1, reads_nc_rep1, tpm_nc_1_rep1, reads_nc_1_rep1, tpm_nc_m_rep1, reads_nc_m_rep1, labels  = run(file_dir+'cp/', 'Vero_mapped_24-1', dict_gene)
tpm_c_rep2, reads_c_rep2, tpm_nc_rep2, reads_nc_rep2, tpm_nc_1_rep2, reads_nc_1_rep2, tpm_nc_m_rep2, reads_nc_m_rep2, labels  = run(file_dir+'cp/', 'Vero_mapped_24-2', dict_gene)
tpm_c_rep3, reads_c_rep3, tpm_nc_rep3, reads_nc_rep3, tpm_nc_1_rep3, reads_nc_1_rep3, tpm_nc_m_rep3, reads_nc_m_rep3, labels  = run(file_dir+'cp/', 'Vero_mapped_24-3', dict_gene)

tpm_c_mean, tpm_c_std = np.mean([tpm_c_rep1, tpm_c_rep2, tpm_c_rep3], axis=0), np.std([tpm_c_rep1, tpm_c_rep2, tpm_c_rep3], axis=0)
reads_c_mean, reads_c_std = np.mean([reads_c_rep1, reads_c_rep2, reads_c_rep3], axis=0), np.std([reads_c_rep1, reads_c_rep2, reads_c_rep3], axis=0)
tpm_nc_mean, tpm_nc_std = np.mean([tpm_nc_rep1, tpm_nc_rep2, tpm_nc_rep3], axis=0), np.std([tpm_nc_rep1, tpm_nc_rep2, tpm_nc_rep3], axis=0)
reads_nc_mean, reads_nc_std = np.mean([reads_nc_rep1, reads_nc_rep2, reads_nc_rep3], axis=0), np.std([reads_nc_rep1, reads_nc_rep2, reads_nc_rep3], axis=0)
tpm_nc_1_mean, tpm_nc_1_std = np.mean([tpm_nc_1_rep1, tpm_nc_1_rep2, tpm_nc_1_rep3], axis=0), np.std([tpm_nc_1_rep1, tpm_nc_1_rep2, tpm_nc_1_rep3], axis=0)
reads_nc_1_mean, reads_nc_1_std = np.mean([reads_nc_1_rep1, reads_nc_1_rep2, reads_nc_1_rep3], axis=0), np.std([reads_nc_1_rep1, reads_nc_1_rep2, reads_nc_1_rep3], axis=0)
tpm_nc_m_mean, tpm_nc_m_std = np.mean([tpm_nc_m_rep1, tpm_nc_m_rep2, tpm_nc_m_rep3], axis=0), np.std([tpm_nc_m_rep1, tpm_nc_m_rep2, tpm_nc_m_rep3], axis=0)
reads_nc_m_mean, reads_nc_m_std = np.mean([reads_nc_m_rep1, reads_nc_m_rep2, reads_nc_m_rep3], axis=0), np.std([reads_nc_m_rep1, reads_nc_m_rep2, reads_nc_m_rep3], axis=0)

n = len(labels)
rt = 90
plt.figure(figsize=(10,4))

plt.subplot(241)
plt.bar(np.arange(n), tpm_c_mean, yerr = tpm_c_std, capsize = 3, color = 'lightblue')
plt.xticks(np.arange(n), labels, rotation = rt)
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)
plt.ylabel('TPM')
plt.title('Canonical')
plt.subplot(245)
plt.bar(np.arange(n), reads_c_mean, yerr = reads_c_std, capsize = 3, color = 'lightblue')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('Reads')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)

plt.subplot(242)
plt.bar(np.arange(n), tpm_nc_mean, yerr = tpm_nc_std, capsize = 3, color = 'grey')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('TPM')
plt.title('none-Canonical')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)
plt.subplot(246)
plt.bar(np.arange(n), reads_nc_mean, yerr = reads_nc_std, capsize = 3, color = 'grey')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('Reads')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)

plt.subplot(243)
plt.bar(np.arange(n), tpm_nc_1_mean, yerr = tpm_nc_1_std, capsize = 3, color = 'pink')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('TPM')
plt.title('once-junction')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)
plt.subplot(247)
plt.bar(np.arange(n), reads_nc_1_mean, yerr = reads_nc_1_std, capsize = 3, color = 'pink')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('Reads')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)

plt.subplot(244)
plt.bar(np.arange(n), tpm_nc_m_mean, yerr = tpm_nc_m_std, capsize = 3, color = 'pink')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('TPM')
plt.title('multiple-junction')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)
plt.subplot(248)
plt.bar(np.arange(n), reads_nc_m_mean, yerr = reads_nc_m_std, capsize = 3, color = 'pink')
plt.xticks(np.arange(n), labels, rotation = rt)
plt.ylabel('Reads')
ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)

tpm_24_nc_m = [24.47, 28.02, 26.47]
tpm_48_nc_m = [36.70, 37.86, 35.05]

tpm_24_nc_1 = [ 430.17, 464.58, 450.75]
tpm_48_nc_1 = [ 267.66, 273.32, 271.22]

sns.set_theme(style="whitegrid")
tpm_24 = tpm_24_nc_m.copy()
tpm_48 = tpm_48_nc_m.copy()
fig = plt.figure(figsize=(6,6))
plt.bar(np.arange(3), tpm_24, color = 'pink', label = '24 hpi')
plt.bar(np.arange(4,7), tpm_48, color = 'lightblue', label = '48 hpi')
plt.legend()
plt.xticks(range(7), ['rep1','rep2','rep3','','rep1','rep2','rep3'])
plt.title(str(ss.ttest_ind(tpm_24, tpm_48)[1]/2))

