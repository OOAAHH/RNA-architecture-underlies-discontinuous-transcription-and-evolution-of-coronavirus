# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:12:36 2022

@author: DELL
"""
from collections import OrderedDict, Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sys import argv

pyname, chrs, file_in, file_name = argv 

# dict_gene = OrderedDict([('leader', (0, 100)),
#  ('UTR_5', (101, 265)),
#  ('ORF1a', (266, 13468)),
#  ('ORF1b', (13468, 21555)),
#  ('S', (21563, 25384)),
#  ('ORF3a', (25393, 26220)),
#  ('E', (26245, 26492)),
#  ('M', (26503, 27191)),
#  ('ORF6', (27202, 27387)),
#  ('ORF7a', (27394, 27759)),
#  ('ORF7b', (27756, 27887)),
#  ('ORF8', (27894, 28259)),
#  ('N', (28274, 29533)),
#  ('ORF9b', (28284, 28577)),
#  ('ORF9c', (28734, 28955)),	
#  ('ORF10', (29558, 29674)),
#  ('UTR_3', (29675, 29903))])

dict_gene = OrderedDict([('leader', (0, 100)),
 ('UTR_5', (101, 539)),
 ('ORF1a', (540, 11414)),
 ('ORF1b', (11414, 19342)),
 ('S', (19324, 22803)),
 ('E', (22797, 23048)),
 ('M', (23041, 23694)),
 ('NS6', (23694, 23978)),
 ('N', (23999, 25027)),
 ('NS7', (24093, 24695)),
 ('New1', (24194, 24659)),
 ('New2', (25040, 25240)),
 ('UTR_3', (25241, 25420))])

#genes = ['ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b', 'ORF9c', 'ORF10']
genes = list(dict_gene.keys())[4:-1]

# chrs = 'MK584552'
# file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/PEDV/nanopore/'
# file_name = 'MK584552.1.raw.splice'
# dict_gene = OrderedDict([('leader', (1, 100)),
#              ('UTR_5', (101, 292)),
#              ('ORF1a', (293, 12601)),
#              ('ORF1b', (12601, 20637)),
#              ('S', (20634, 24791)),
#              ('ORF3', (24791, 25465)),
#              ('E', (25446, 25676)),
#              ('M', (25684, 26364)),
#              ('N', (26376, 27701)),
#              ('N1', (26836, 27024)),
#              ('New', (27782, 27919)),
#              ('UTR_3', (27920, 28044))])
# genes = list(dict_gene.keys())[4:-1]

def is_in(dict_data, in_data):
    out_data = []
    for key in dict_data:
        values = dict_data[key]
        if values[0] <=in_data and values[1] >= in_data:
            out_data.append(key)
    return out_data

def is_between(dict_data, in_data):
    out_data = []
    for i,key in enumerate(dict_data):
        values = dict_data[key]
        
        if i == 0:
            if values[0] <=in_data[0] and values[1] >= in_data[0]:
                out_data.append(key)
        
        else:
            if values[0] >=in_data[0] and values[1] <= in_data[1]:
                out_data.append(key)
                
    return out_data

def is_leader(dict_data, in_data):
    out_data = []
    for i,key in enumerate(dict_data):
        values = dict_data[key]
        
        if i == 0:
            if values[0] <=in_data[0] and values[1] >= in_data[0]:
                out_data.append(key)
        
        else:
            if values[0] >=in_data[0] and values[1] <= in_data[1]:
                out_data.append(key)
                
    return out_data

def junction(dict_data, in_data):
    out_data = []
    for key in dict_data:
        values = dict_data[key]
        if values[0] <=in_data[0] and values[1] >= in_data[0]:
            out_data.append(key)
        
    if len(out_data) == 0:    
        out_data.append('gap')
        
    for key in dict_data:
        values = dict_data[key]        
        if values[0] <=in_data[1] and values[1] >= in_data[1]:
            out_data.append(key)
            
    if len(out_data) == 1:    
            out_data.append('gap')
            
    return '-'.join(out_data)

def delete_list(list_in, list_delet):
    for x in list_delet:
        if x in list_in:
            list_in.remove(x)
    return list_in
    

files_in = open(file_in)
files_out2 = open(file_name + '-2.bed', 'w')
files_out3 = open(file_name + '-3.bed', 'w')
transcpt_out = {}
transcpt_out_none = {}
reads_count = 0
max_junc = 0
for line in files_in:
    lines = line.split('\t')  
    if lines[0] >= 125:
        continue
    reads_count+=1
    lens = len(lines)
    iters = int(lens/2)
    names = []
    trscpt = []
    trscpt_junc = []
    if lens % 2 != 0:
        print('die out')
        break
    for i in range(iters):
        start = int(lines[i*2])
        end = int(lines[i*2+1])
        if i == (iters-1) :
            end = dict_gene['UTR_3'][1] 
        gene_name = is_between(dict_gene, [start, end])
        trscpt.extend(gene_name)
        trscpt_junc.append(junction(dict_gene,[start, end]))
        names.extend( [str(start), str(end)]) 
        
    trscpt = delete_list(trscpt, ['leader', 'UTR_5'])
        
    if len(trscpt) == 0 :
        continue
    juncs = iters - 1
    if juncs > max_junc :
        max_junc = juncs
        
    true_gen = trscpt[0]
    keys = (true_gen, juncs)
    values = '->'.join(trscpt_junc)
    
    if trscpt_junc[0].split('-')[1] == 'leader':
        num = len(trscpt_junc)
        if num == 2 :
            if transcpt_out.__contains__(keys):
                transcpt_out[keys].append(values)
            else:
                transcpt_out[keys] = []
                transcpt_out[keys].append(values)
        
        if num > 2:
            if transcpt_out_none.__contains__(keys):
                transcpt_out_none[keys].append(values)
            else:
                transcpt_out_none[keys] = []
                transcpt_out_none[keys].append(values)
    else:
        if transcpt_out_none.__contains__(keys):
            transcpt_out_none[keys].append(values)
        else:
            transcpt_out_none[keys] = []
            transcpt_out_none[keys].append(values)
                
    # if true_gen == 'ORF9b':
    if len(names) >= 5:
        files_out3.writelines('\t'.join(names)+'\t'+values+'\n')
    else :
        files_out2.writelines('\t'.join(names)+'\t'+values+'\n')
    # if len(bed) >= 7:
    #     files_out.writelines('\t'.join(bed)+'\n')
files_out2.close()
files_out3.close()
files_in.close()
#print(Counter(transcpt_out[('ORF10', 1)]))

transcpt_out_count = {}
transcpt_out_none_count = {}
for gene in genes:
    for i in range(max_junc):
        keys = (gene, i+1)
        if transcpt_out.__contains__(keys):
              reads_count =  sum(Counter(transcpt_out[keys]).values())
              if transcpt_out_count.__contains__(gene):
                  transcpt_out_count[gene][0] += reads_count
              else :
                  transcpt_out_count[gene] = [0,0]
                  transcpt_out_count[gene][0] += reads_count

        if transcpt_out_none.__contains__(keys) and i == 0:
              reads_count =  sum(Counter(transcpt_out_none[keys]).values())
              if transcpt_out_none_count.__contains__(gene):
                  transcpt_out_none_count[gene][0] += reads_count

              else :
                  transcpt_out_none_count[gene] = [0,0,0,0]
                  transcpt_out_none_count[gene][0] += reads_count
        
        if transcpt_out_none.__contains__(keys) and i != 0:
              reads_count =  sum(Counter(transcpt_out_none[keys]).values())
              if transcpt_out_none_count.__contains__(gene):
                  transcpt_out_none_count[gene][2] += reads_count

              else :
                  transcpt_out_none_count[gene] = [0,0,0,0]
                  transcpt_out_none_count[gene][2] += reads_count

genes = [key for key in transcpt_out_count.keys()]                  
genes_none = [key for key in transcpt_out_none_count.keys()] 
gene_lens = [dict_gene[g][1] - dict_gene[g][0] for g in genes]
gene_none_lens = [dict_gene[g][1] - dict_gene[g][0] for g in genes_none]
  
none_R_L = sum([transcpt_out_none_count[g][0]/gene_none_lens[i] for i,g in enumerate(genes_none)])
R_L = sum([transcpt_out_count[g][0]/gene_lens[i] for i,g in enumerate(genes)]) 

sum_R_L = none_R_L + R_L
for i,gene in enumerate(genes):
    transcpt_out_count[gene][1] = transcpt_out_count[gene][0]  / sum_R_L / gene_lens[i] * 1000

for i,gene in enumerate(genes_none):
    transcpt_out_none_count[gene][1] = transcpt_out_none_count[gene][0] / sum_R_L / gene_none_lens[i] * 1000
    transcpt_out_none_count[gene][3] = transcpt_out_none_count[gene][2] / sum_R_L / gene_none_lens[i] * 1000

print(sum_R_L)
print(transcpt_out_none_count)
print(transcpt_out_count)

num = len(genes)
num_none = len(genes_none)
rt = 90
fig = plt.figure(figsize=(12, 8.5))

plt.subplot(241)
plt.bar(np.arange(num),[x[1] for x in transcpt_out_count.values()], color = 'lightblue')
plt.xticks(np.arange(num),[x for x in transcpt_out_count.keys()], rotation = rt)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.ylabel('TPM')
plt.title('Canonical')
plt.subplot(245)
plt.bar(np.arange(num),[x[0] for x in transcpt_out_count.values()], color = 'lightblue')
plt.xticks(np.arange(num),[x for x in transcpt_out_count.keys()], rotation = rt)
plt.ylabel('Reads')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.subplot(242)
plt.bar(np.arange(num_none),[x[1] + x[3] for x in transcpt_out_none_count.values()], color = 'grey')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('TPM')
plt.title('none-Canonical')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplot(246)
plt.bar(np.arange(num_none),[x[0] + x[2] for x in transcpt_out_none_count.values()], color = 'grey')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('Reads')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.subplot(243)
plt.bar(np.arange(num_none),[x[1] for x in transcpt_out_none_count.values()], color = 'pink')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('TPM')
plt.title('one-juctiong')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplot(247)
plt.bar(np.arange(num_none),[x[0] for x in transcpt_out_none_count.values()], color = 'pink')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('Reads')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.subplot(244)
plt.bar(np.arange(num_none),[x[3] for x in transcpt_out_none_count.values()], color = 'pink')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('TPM')
plt.title('multiple-juction')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplot(248)
plt.bar(np.arange(num_none),[x[2] for x in transcpt_out_none_count.values()], color = 'pink')
plt.xticks(np.arange(num_none),[x for x in transcpt_out_none_count.keys()], rotation = rt)
#plt.ylabel('Reads')
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

pp = PdfPages(file_name + '-transcription.pdf')
pp.savefig(fig)
pp.close()
plt.close()
##           
# res = 10 
# data = np.loadtxt(file_dir + file_name + '-2.bed', usecols=[0,1,2,3,], dtype=[('S1','<i4'),('E1','<i4'),('S2','<i4'),('E2','<i4')],)
# ma_junc_2 = np.asarray(list(zip(data['E1'] / 10,data['S2'] / 10,np.ones(data.shape[0]))), dtype=[('S','i4'),('E','i4'),('value','f8')])    
# ma_junc_2 = get_matrix_npz(ma_junc_2, ma_junc_2)
# plt.figure()
# plt.imshow(ma_junc_2[0:500, 2500:], cmap= plt.get_cmap('Reds'), vmax=10)
