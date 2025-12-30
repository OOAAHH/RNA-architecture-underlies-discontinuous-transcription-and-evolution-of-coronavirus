# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 10:58:17 2023

@author: Administrator
"""

import matplotlib.pyplot as plt
import numpy as np
from SARS2_sgRNA import rna_seq_24_all, nanopore_24_all, nanopore_48_all, nanopore_kim_24
from genome_382_short import dict_gene_382, dict_gene_346, dict_gene_220, dict_gene_188, dict_gene_142, dict_gene_133, dict_gene_90, dict_gene_62, dict_gene, dict_gene_2
from genome_382_short import dict_gene_186, dict_gene_171, dict_gene_207, dict_gene_154, dict_gene_116, dict_gene_92, dict_gene_78, dict_gene_50, dict_gene_78
from Coronavirus_ct_shannon_data import strc_ric_sars2_virion

res = 1000
color = ['lightgreen','lightgreen', 'pink', 'lightblue', 'pink', 'pink', 'lightblue', 'lightblue', 'lightblue', 'lightblue', 'pink', 'lightblue', 'lightblue', 'lightblue']
color_345 = ['lightblue', 'lightblue', 'lightblue', 'lightblue', 'lightblue', 'pink', 'lightblue', 'lightblue', 'lightblue']

recom_strc = np.zeros(strc_ric_sars2_virion.shape[0], dtype=[('S', '<i4'), ('E', '<i4'), ('value', '<f8')])
recom_strc['S'] = strc_ric_sars2_virion['1']
recom_strc['E'] = strc_ric_sars2_virion['5']
recom_strc['value'] = strc_ric_sars2_virion['5']
recom_strc_hp = recom_strc[recom_strc['value']>0]
recom_strc_hp['value'] = 100

recom_data = rna_seq_24_all #recom_strc_hp #
distance = recom_data['E'] - recom_data['S']
recom_data_long = recom_data[distance >= 100]
recom_data_short = recom_data[distance <= 1000]
recom_data_highq = recom_data_short[ recom_data_short['value'] > 20]
start = 27790
end = 28229 + 50 #end = 29550 #
target_seq_1 = recom_data_long[ (recom_data_long['S'] > start) & (recom_data_long['E'] < end)]
start = 27472 - 30
end = 27760
target_seq_2 = recom_data_long[ (recom_data_long['S'] > start) & (recom_data_long['E'] < end)]
color_reads = plt.get_cmap('Greys')( np.linspace( 0.4, 0.8, int( np.log2( max(target_seq_1['value'].max(), target_seq_2['value'].max() ) ) ) + 1 ) )
color_reads_hq = plt.get_cmap('Greys')(np.linspace(0.1, 0.9, int(np.log2(recom_data_highq['value'].max())) + 1))

genome_id = ['382', '346', '188', '142','133','90', '62', '220',
             '186', '171', '207', '154', '116', '92', '78', '50']
list_id = [ '346','220', '188','62','186', '171', '207', '154', '116', '92', '78', '50']

fig = plt.figure(figsize=(14,8))
# plt.subplot(211)
i = 6
for key in dict_gene_2:
    s,e = dict_gene_2[key]
    s = s / res
    e = e / res
    if key.split('_')[0] == 'gap':
        plt.plot([s,e],[0,0], color = 'black', alpha = 0.5)
        
    else:
        plt.fill_between([s,e],[0.05,0.05], -0.05, color = color[i], alpha = 0.5)
        i += 1

for id_seq in ['1','2']:
    target_seq = locals()['target_seq_' + id_seq]
    for seq in target_seq[:]:
        s = seq[0] / res
        e = seq[1] / res
        r = (e - s) / 2
        a,b = ((s+e)/2, 0.05)
        value = int(np.log2(seq[2]))
        x = np.arange(a - r, a + r, 0.0001 )
        y = b + np.sqrt(r ** 2 - (x - a) ** 2)
        plt.plot(x , y, color = color_reads[value], lw = 0.15, alpha = 0.5)

up = -0.07
for id_g in genome_id:
    genome = locals()['dict_gene_' + id_g]
    i = 6
    if id_g in list_id:
         color_genome = color_345
         i = 0
    else :
        color_genome = color
    for key in genome:
        s,e = genome[key]
        s = s / res
        e = e / res    
        down = up - 0.1
        mid = (up + down) / 2
        if key.split('_')[0] == 'gap':
            plt.plot([s,e],[mid, mid], color = 'black', alpha = 0.5)
            
        elif key.split('_')[0] == 'delta':
             plt.plot([s,e],[mid, mid], color = 'red', alpha = 0.5)    
             
        else:
            plt.fill_between([s,e],[up, up], down, color = color_genome[i], alpha = 0.5)
            i += 1
    up = down - 0.02

second = recom_strc_hp[(recom_strc_hp['S'] > 27250) & (recom_strc_hp['E'] < 29500)]
for seq in second[second['S'] < second['E']]:
    s = seq[0] / res
    e = seq[1] / res
    r = (e - s) / 2
    a,b = ((s+e)/2, down)
    value = int(np.log2(seq[2]))
    x = np.arange(a - r, a + r, 0.0001 )
    y = b + np.sqrt(r ** 2 - (x - a) ** 2)
    plt.plot(x , -y + down*2, color = color_reads_hq[value], lw = 0.15)   

plt.xlim(27.25, 29.5)
plt.ylim(-3.125, 0.35)
plt.yticks([])
ax = plt.gca()
ax.spines[['left', 'right', 'top']].set_visible(False)


