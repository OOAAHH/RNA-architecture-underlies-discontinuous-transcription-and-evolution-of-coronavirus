# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 11:01:34 2023

@author: Administrator
"""
from SARS2_deletion import sars2_long_del, sars2_rand_del
from PEDV_deletion import pedv_long_del, pedv_rand_del
from SARS2_recom_data import rna_seq_24_all
from Coronavirus_recom_data import rna_pedv_cell_all
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss

sras2_recomb = rna_seq_24_all[rna_seq_24_all['value']>1].copy()
sars2_true = []
sars2_ran = []
for i in np.arange(len(sars2_long_del)):
    true_del = sars2_long_del[i]
    ran_del = sars2_rand_del[i]
    distance_del = abs(sras2_recomb['S'] - true_del[0]) + abs(sras2_recomb['E'] - true_del[1])
    distance_ran = abs(sras2_recomb['S'] - ran_del[0]) + abs(sras2_recomb['E'] - ran_del[1])
    if distance_del.min() > 100:
        print(true_del)
    sars2_true.append(distance_del.min())
    sars2_ran.append(distance_ran.min())
    
PEDV_recomb = rna_pedv_cell_all[rna_pedv_cell_all['value'] > 1].copy()
pedv_true = []
pedv_ran = []
for i in np.arange(len(pedv_long_del)):
    true_del = pedv_long_del[i]
    ran_del = pedv_rand_del[i]
    distance_del = abs(PEDV_recomb['S'] - true_del[0]) + abs(PEDV_recomb['E'] - true_del[1])
    distance_ran = abs(PEDV_recomb['S'] - ran_del[0]) + abs(PEDV_recomb['E'] - ran_del[1])
    if distance_del.min() > 100:
        print(true_del)
    if distance_del.min() < 5:
        print(true_del, PEDV_recomb[distance_del<5])
    pedv_true.append(distance_del.min())
    pedv_ran.append(distance_ran.min())    

df_sars2 = pd.DataFrame({'True':sars2_true,'Random':sars2_ran})
df_pedv = pd.DataFrame({'True':pedv_true,'Random':pedv_ran})
fig = plt.figure(figsize=(12,6))
plt.subplot(223)
sns.boxplot(data=df_sars2)
plt.text(0.15, 100, 'P-value:%.2e'%ss.ranksums(sars2_true, sars2_ran, alternative='less')[1])
#sns.stripplot(data=df_sars2)

plt.subplot(221)
sns.histplot(data=df_sars2, bins=20, binrange=(1,100), alpha = 0.5, kde = True)
plt.title('SARS2')

plt.subplot(224)
sns.boxplot(data=df_pedv)
plt.text(0.15, 100, 'P-value:%.2e'%ss.ranksums(pedv_true, pedv_ran, alternative='less')[1])

plt.subplot(222)
sns.histplot(data=df_pedv, bins=20, binrange=(1,100), alpha = 0.5, kde = True)
plt.title('PEDV')


