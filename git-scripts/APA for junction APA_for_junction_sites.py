# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:11:19 2023

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt
from ric_fun import  apa, classify_rna_data, get_oe
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from SARS2_sgRNA import rna_seq_24_all, nanopore_24_all, nanopore_48_all, nanopore_kim_24
from Coronavirus_hic_data import dis_sars_rna_cell
from Coronavirus_hic_data import ma_sars_com_cell, oe_sars_com_cell, dis_sars_com_cell, ma_sars_splash_cell, oe_sars_splash_cell, dis_sars_splash_cell
from Coronavirus_hic_data import ma_sars_ric_virion_kr, oe_sars_ric_virion_kr, dis_sars_ric_virion
import scipy.stats as ss
oe_com, dis_com = get_oe(ma_sars_com_cell)

pp = PdfPages('Figure_2-SARS2-APA.pdf')
res = 5
###距离衰减曲线
#sns.set_theme(style="darkgrid")
sns.set_theme(style="ticks")
fig = plt.figure(figsize=(8,4))
s = 300
e = 3000
s = int(s/ res)
e = int(e/ res) 
lw = 2
#plt.plot(dis_sars_com_cell[s:e], label = 'SPLASH')
plt.plot(np.arange(s,e), dis_sars_rna_cell[s:e], lw = lw, color = 'black', label = 'RNA-seq')
plt.plot(np.arange(s,e), dis_sars_ric_virion[s:e], lw = lw, color = 'green', label = 'vRIC-seq')
plt.plot(np.arange(s,e), dis_sars_splash_cell[s:e], lw = lw, color = 'grey', label = 'SPLASH')
plt.plot(np.arange(s,e), dis_com[s:e], lw = lw, color = 'lightblue', label = 'COMRADES')
plt.xticks( np.arange(s,e,100), np.arange(s,e,100)*res)
plt.ylabel('Number of chimeric reads')
plt.xlabel('Distance between junction sites (bp)', labelpad = 0.1)
plt.legend()
ax = plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.annotate('P1', xy=(160, 0.25), xycoords='data',
            xytext=(120, 1.6), textcoords='data',
            va='top', ha='left',
            arrowprops=dict(facecolor='black', shrink=0.01))
ax.annotate('P2', xy=(277, 0.375), xycoords='data',
            xytext=(260, 2.6), textcoords='data',
            va='top', ha='left',
            arrowprops=dict(facecolor='black', shrink=0.01))
ax.annotate('P3', xy=(420, 0.175), xycoords='data',
            xytext=(450, 1.45), textcoords='data',
            va='top', ha='left',
            arrowprops=dict(facecolor='black', shrink=0.01))

pp.savefig(fig)
plt.close()

##center data
box_ric_short = []
box_ric_mid = []
box_ric_long = []
box_ric_can = []
box_com_short = []
box_com_mid = []
box_com_long = []
box_com_can = []
box_sp_short = []
box_sp_mid = []
box_sp_long = []
box_sp_can = []
recom_data = rna_seq_24_all.copy() #nanopore_kim_24.copy() #
rna_c, rna_nc, hist_data = classify_rna_data(recom_data, 100, 10, 100, [100, 2000])
### APA
res = 5
distance = rna_nc['E'] - rna_nc['S']
fig = plt.figure(figsize=(16,14))
ma_in = oe_sars_ric_virion_kr.copy()
fig.subplots_adjust(wspace=0.1, hspace=0.2)
sns.set_theme(style="whitegrid")
plt.subplot(321)
plt.plot(hist_data, color = 'grey')
#plt.fill_between([25,200], [9800,9800], -300, alpha = 0.25, color = 'grey')
plt.plot([25,25], [-300,9800], linestyle='dashed', color = 'pink')
plt.plot([200,200], [-300,9800], linestyle='dashed', color = 'pink')
plt.xticks(np.arange(0,300,50), np.arange(0,300,50)*100)
plt.xlabel('Distance between junction sites (bp)')
plt.ylabel('Number of chimeric reads')

plt.subplot(322)
label = ['< 2500 bp','2500-20000 bp','> 20000 bp']
colors = plt.get_cmap('Reds')(np.linspace(0.2, 0.7, 3))
data_in = [sum(hist_data[:25]),sum(hist_data[25:200]),sum(hist_data[200:])]
plt.pie(data_in, colors=colors, radius=1.35, labeldistance=1.15, autopct='%.2f%%', center=(4, 4),
       wedgeprops={"linewidth": 1, "edgecolor": "white"}, labels = label )

sns.set_theme(style="ticks")
plt.subplot(6,4,9)
ma, count, data = apa(ma_in, rna_nc[(distance  > 20000)], 0, 30000, 10, res)
box_ric_long.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('distance > 20000 bp')
plt.ylabel('vRic-seq')

plt.subplot(6,4,10)
ma, count, data = apa(ma_in, rna_nc[(distance  < 20000) & (distance  > 2500)], 0, 30000, 10, res)
box_ric_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 35, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('2500bp < distance < 20000 bp')

plt.subplot(6,4,11)
ma, count, data = apa(ma_in, rna_nc[(distance  < 2500)], 0, 30000, 10, res, 'stem')
box_ric_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 75, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('distance < 2500 bp')

plt.subplot(6,4,12)
ma, count, data = apa(ma_in, rna_c, 0, 30000, 10, res, 'diag')
box_ric_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 35, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('canonical subgenome')

ma_in = oe_sars_splash_cell.copy()
plt.subplot(6,4,13)
ma, count, data = apa(ma_in, rna_nc[(distance  > 20000)], 0, 30000, 10, res)
box_sp_long.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.ylabel('SPLASH')

plt.subplot(6,4,14)
ma, count, data = apa(ma_in, rna_nc[(distance  < 20000) & (distance  > 2500)], 0, 30000, 10, res)
box_sp_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(6,4,15)
ma, count, data = apa(ma_in, rna_nc[(distance  < 2500)], 0, 30000, 10, res, 'stem')
box_sp_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(6,4,16)
ma, count, data = apa(ma_in, rna_c, 0, 30000, 10, res, 'diag')
box_sp_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

ma_in = oe_com.copy()
plt.subplot(6,4,17)
ma, count, data = apa(ma_in, rna_nc[(distance  > 20000)], 0, 30000, 10, res)
box_com_long.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.ylabel('COMRADES')

plt.subplot(6,4,18)
ma, count, data = apa(ma_in, rna_nc[(distance  < 20000) & (distance  > 2500)], 0, 30000, 10, res)
box_com_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(6,4,19)
ma, count, data = apa(ma_in, rna_nc[(distance  < 2500)], 0, 30000, 10, res, 'stem')
box_com_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(6,4,20)
ma, count, data = apa(ma_in, rna_c, 0, 30000, 10, res, 'diag')
box_com_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
pp.savefig(fig)
plt.close()


fig = plt.figure(figsize=(14,8))
lens = rna_nc[distance  > 20000].shape[0]
rs = np.random.random_sample(lens * 5) * 9903 
re = rs + np.random.random_sample(lens * 5) * 9903 + 20000
mask = re < 29903
rand = rna_nc[distance  > 20000].copy()
rand['S'] = rs[mask][:lens]
rand['E'] = re[mask][:lens]

ma_in = oe_sars_ric_virion_kr.copy()
ma, count, data = apa(ma_in, rand, 50, 30000, 10, res)
box_ric_long.append(data)
plt.subplot(3,4,1)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('distance > 20000 bp')
plt.ylabel('vRic-seq')

ma_in = (oe_sars_splash_cell.copy())
plt.subplot(3,4,5)
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res)
box_sp_long.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.ylabel('SPLASH')

ma_in = (oe_com.copy())
plt.subplot(3,4,9)
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res)
box_com_long.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.ylabel('COMRADES')

rs = np.random.random_sample(rna_nc[(distance  < 20000) & (distance  > 2500)].shape[0]) * 10000 
re = rs + np.random.random_sample(rna_nc[(distance  < 20000) & (distance  > 2500)].shape[0]) * 17500 + 2403
rand = rna_nc[(distance  < 20000) & (distance  > 2500)].copy()
rand['S'] = rs
rand['E'] = re

plt.subplot(3,4,2)
ma_in = oe_sars_ric_virion_kr.copy()
ma, count, data = apa(ma_in, rand, 50, 30000, 10, res)
box_ric_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 35, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('2500bp < distance < 20000 bp')

ma_in = oe_sars_splash_cell.copy()
plt.subplot(3,4,6)
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res)
box_sp_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

ma_in = oe_com.copy()
plt.subplot(3,4,10)
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res)
box_com_mid.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

rs = np.random.random_sample(rna_nc[(distance  < 2500)].shape[0]) * 27403
re = rs + np.random.random_sample(rna_nc[(distance  < 2500)].shape[0]) * 2500
rand = rna_nc[(distance  < 2500)].copy()
rand['S'] = rs
rand['E'] = re

plt.subplot(3,4,3)
ma_in = oe_sars_ric_virion_kr.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'stem')
box_ric_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 75, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('distance < 2500 bp')

plt.subplot(3,4,7)
ma_in = oe_sars_splash_cell.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'stem')
box_sp_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(3,4,11)
ma_in = oe_com.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'stem')
box_com_short.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 20, vmin = 0.5)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

rs =  np.random.random_sample(rna_c.shape[0]) * 100
re = rs + np.random.random_sample(rna_c.shape[0]) * 29803
rand = rna_c.copy()
rand['S'] = rs
rand['E'] = re

plt.subplot(3,4,4)
ma_in = oe_sars_ric_virion_kr.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'diag')
box_ric_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
plt.title('canonical subgenome')

plt.subplot(3,4,8)
ma_in = oe_sars_splash_cell.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'diag')
box_sp_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)

plt.subplot(3,4,12)
ma_in = oe_com.copy()
ma, count, data = apa(ma_in, rand, 0, 30000, 10, res, 'diag')
box_com_can.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 10, vmin = 0)
#plt.colorbar(location = 'top', shrink = 0.5)
plt.xticks([0,10,20],["5'",'junction sites',"3'"])
plt.yticks([0,10,20],["5'",'junction\n sites',"3'"], rotation = 0)
pp.savefig(fig)
plt.close()

fig = plt.figure(figsize=(12,6))
sns.set_theme(style="whitegrid")
plt.subplot(341)
sns.boxplot(data=box_ric_long, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_ric_long[0], box_ric_long[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,150])
plt.xticks([0,1],['True','Random'])
plt.ylabel('Interaction Strength')
plt.title('distance > 20000 bp')

plt.subplot(342)
sns.boxplot(data=box_ric_mid, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_ric_mid[0], box_ric_mid[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,150])
plt.xticks([0,1],['True','Random'])
plt.title('2500bp < distance < 20000 bp')

plt.subplot(343)
sns.boxplot(data=box_ric_short, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_ric_short[0], box_ric_short[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,350])
plt.xticks([0,1],['True','Random'])
plt.title('distance < 2500 bp')

plt.subplot(344)
sns.boxplot(data=box_ric_can, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_ric_can[0], box_ric_can[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,110])
plt.xticks([0,1],['True','Random'])
plt.title('Canonical subgenome')

plt.subplot(345)
sns.boxplot(data=box_sp_long, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_sp_long[0], box_sp_long[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,100])
plt.xticks([0,1],['True','Random'])
plt.ylabel('Interaction Strength')

plt.subplot(346)
sns.boxplot(data=box_sp_mid, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_sp_mid[0], box_sp_mid[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,110])
plt.xticks([0,1],['True','Random'])

plt.subplot(347)
sns.boxplot(data=box_sp_short, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_sp_short[0], box_sp_short[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,350])
plt.xticks([0,1],['True','Random'])

plt.subplot(348)
sns.boxplot(data=box_sp_can, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_sp_can[0], box_sp_can[1],alternative='greater')[1]
plt.text(0.125, 10, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,50])
plt.xticks([0,1],['True','Random'])

plt.subplot(349)
sns.boxplot(data=box_com_long, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_com_long[0], box_com_long[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,100])
plt.xticks([0,1],['True','Random'])
plt.ylabel('Interaction Strength')

plt.subplot(3,4,10)
sns.boxplot(data=box_com_mid, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_com_mid[0], box_com_mid[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,110])
plt.xticks([0,1],['True','Random'])

plt.subplot(3,4,11)
sns.boxplot(data=box_com_short, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_com_short[0], box_com_short[1],alternative='greater')[1]
plt.text(0.125, 50, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,350])
plt.xticks([0,1],['True','Random'])

plt.subplot(3,4,12)
sns.boxplot(data=box_com_can, color='grey') #,showfliers = False
p_value = ss.ttest_ind(box_com_can[0], box_com_can[1],alternative='greater')[1]
plt.text(0.125, 10, 'P-value=%.2e'%p_value) 
#plt.ylim([-1.5,50])
plt.xticks([0,1],['True','Random'])

pp.savefig(fig)
plt.close()
pp.close()


