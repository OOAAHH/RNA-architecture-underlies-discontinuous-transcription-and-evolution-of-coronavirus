# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 11:19:13 2023

@author: Administrator
"""
import matplotlib.pyplot as plt
from collections import Counter
from ric_fun import apa
import numpy as np
from Coronavirus_hic_data import oe_sars_ric_virion
import scipy.stats as ss
import seaborn as sns

dir_pre = 'G:/'
file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/data/mafft/SARS-COV2/'
virus_bg = list(next(open((file_dir + '1110-mafft-done-NC_045512.fasta'))))

res = 5
ma_ric = oe_sars_ric_virion.copy()


dict_loci = {}
j = 0
for i,x in enumerate(virus_bg):
    if x != '-':
        j += 1
    dict_loci[i] = j

def sites_to_range(input_list):
    if len(input_list) == 0:
        return input_list
    else :
        out_list = []
        start = input_list[0]
        tem = start
        for x in input_list[1:]:
            if (x - tem) == 1:
                tem = x
                continue
            else:
                out_list.append([start, tem, tem-start+1])
                start = x
                tem = x
        out_list.append((start, x, x - start + 1))
        return out_list
    
    
file_in =  open(file_dir + '1110-mafft-done.fasta')
dict_sars2 = {}
for line in file_in:
    lines = list(line)
    if lines[0] == '>':
        key = ''.join(lines[1:-1]).split(' ')[0]
        dict_sars2[key] = [[],[],[],[]]
        if key == 'MN692770.1':
            continue
    else:
        tem_loci = 0
        for i,x in enumerate(lines[:-1]):
            loci = dict_loci[i]
            if loci == tem_loci and x == '-':
                dict_sars2[key][0].append(loci) 
            elif loci != tem_loci and x == '-':
                dict_sars2[key][1].append(loci) 
            elif loci == tem_loci and x != '-':
                dict_sars2[key][2].append(loci) 
            elif loci != tem_loci and x != '-':
                dict_sars2[key][3].append(loci) 
            tem_loci = loci
            
file_in.close()

out_del = []
out_ins = []
for key in dict_sars2:
    if len(dict_sars2[key][1]) > 100:
        out_del.append(key)
    if len(dict_sars2[key][2]) > 100:
        out_ins.append(key)
        
###gap/insertion to name
dict_sars2_name_del = {}
dict_sars2_name_ins = {}

for key in dict_sars2:
    if len(dict_sars2[key][1]) > 1:
        genome_gap = sites_to_range(dict_sars2[key][1])  
        genome_ins = Counter(dict_sars2[key][2])
        for gap in genome_gap:
            if tuple([gap[0], gap[1]]) in dict_sars2_name_del.keys():
                dict_sars2_name_del[tuple([gap[0], gap[1]])].append(key)
            else:
                dict_sars2_name_del[tuple([gap[0], gap[1]])] = []
                dict_sars2_name_del[tuple([gap[0], gap[1]])].append(key)
        if len(genome_ins) > 0:
            for ins in genome_ins:
                ins_key = tuple([ins, genome_ins[ins]])
                if ins_key in dict_sars2_name_ins.keys():
                    dict_sars2_name_ins[ins_key].append(key)
                else:
                    dict_sars2_name_ins[ins_key] = []
                    dict_sars2_name_ins[ins_key].append(key)
                
                

file_out = open(file_dir+'SARS2-genome-del-1110','w')        
for key in sorted(dict_sars2_name_del):
    str_key = '\t'.join([str(key[0]), str(key[1]), str(key[1]-key[0]+1)]) 
    str_value = '\t'.join(dict_sars2_name_del[key])
    strs = '\t'.join([str_key, str_value])
    file_out.writelines(strs+'\n')
file_out.close()

file_out = open(file_dir+'SARS2-genome-ins-1110','w')        
for key in sorted(dict_sars2_name_ins):
    str_key = '\t'.join([str(key[0]), str(key[1])]) 
    str_value = '\t'.join(dict_sars2_name_ins[key])
    strs = '\t'.join([str_key, str_value])
    file_out.writelines(strs+'\n')
file_out.close()

file_out = open(file_dir+'SARS2-genome-del-1110.bedpe','w')        
for key in sorted(dict_sars2_name_del):
    if key[0] > 75 and  key[1] < 29900 and (key[1] - key[0])>=3:
        strs= '\t'.join(['HS1', str(key[0]), str(key[1]),'HS1', str(key[0]), str(key[1])]) 
        file_out.writelines(strs+'\n')
file_out.close()

file_out = open(file_dir+'SARS2-genome-ins-1110.bedpe','w')        
for key in sorted(dict_sars2_name_ins):
    if key[1] >= 3:
        if key[0] < 29780 and key[0] > 200:
            strs= '\t'.join(['HS1', str(key[0]), str(key[0]),'HS1', str(key[0]), str(key[0]+key[1])]) 
            file_out.writelines(strs+'\n')
file_out.close()

##long deletion
list_del = []
for key in dict_sars2:
    if len(dict_sars2[key][1]) > 1:
        list_del.extend(sites_to_range(dict_sars2[key][1]))      
list_del = [tuple(x) for x in list_del]   
del_count = Counter(list_del) 

long_del = []   
for key in del_count.keys():
    if del_count[key] >= 1 and key[2]>30:
        if key[0] < 75 :
            continue
        if key[1] > 29900:
            continue   
        long_del.append(key)

for key in dict_sars2:
    if len(dict_sars2[key][1]) > 100:
        print(key, sites_to_range(dict_sars2[key][1]))    

##random
sars2_rand_del = np.zeros(len(long_del), dtype=[('S', '<i4'), ('E', '<i4'), ('value', '<f8')])
sars2_rand_del['S'] = np.random.random_sample(sars2_rand_del.shape[0]) * 27400
sars2_rand_del['E'] = sars2_rand_del['S'] + np.array([x[2] for x in long_del])        
for i in range(10):     
    rand_i = np.zeros(len(long_del), dtype=[('S', '<i4'), ('E', '<i4'), ('value', '<f8')])
    rand_i['S'] = np.random.random_sample(rand_i.shape[0]) * 27400
    rand_i['E'] = rand_i['S'] + np.array([x[2] for x in long_del])
    sars2_rand_del = np.hstack([sars2_rand_del, rand_i])

sars2_rand_del = sars2_rand_del[:1000]

##APA
sars2_long_del = long_del.copy()   
sns.set_theme(style="ticks")    
plt.figure(figsize=(10,5))
box_data = []
ma_len = 10
plt.subplot(141)
ma, count, data = apa(ma_ric, sars2_long_del, 0, 30000, ma_len, res, 'stem')
box_data.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 5, vmin = 0)
plt.xticks([0,ma_len,2*ma_len],["5'",'juction sites',"3'"])
plt.yticks([0,ma_len,2*ma_len],["5'",'juction\n sites',"3'"], rotation = 0)
plt.title('Genomic deletion')

plt.subplot(142)
ma, count, data = apa(ma_ric, sars2_rand_del, 0, 30000, ma_len, res, 'stem')
box_data.append(data)
plt.imshow(ma / count, cmap = plt.get_cmap('Reds'), vmax = 5, vmin = 0)
plt.xticks([0,ma_len,2*ma_len],["5'",'juction sites',"3'"])
plt.yticks([0,ma_len,2*ma_len],["5'",'juction\n sites',"3'"], rotation = 0)
plt.title('Random')

sns.set_theme(style="whitegrid")
plt.subplot(122)
sns.boxplot(data=box_data, color='grey') #,showfliers = False
p_value = ss.mannwhitneyu(box_data[0], box_data[1],alternative='greater')[1]
plt.text(0.125, 2, 'P-value=%.4e'%p_value) 
plt.ylim([-0.5,15])
plt.xticks([0,1],['Genomic deletion','Random'])
plt.ylabel('Interaction Strength')

##insertion           
list_ins = []
for key in dict_sars2:
    if len(dict_sars2[key][2]) > 1:
        count_ins = Counter(dict_sars2[key][2])
        list_ins.extend([[x, count_ins[x]] for x in count_ins.keys()])      
list_ins = [tuple(x) for x in list_ins]   

ins_count = Counter(list_ins) 
true_ins = []   
for key in ins_count.keys():
    if ins_count[key] >= 1 :
        if key[0] < 205 :
            continue
        if key[0] > 29780:
            continue      
        if key[1] < 3:
            continue
        true_ins.append(key)

sars_ins = true_ins.copy()

