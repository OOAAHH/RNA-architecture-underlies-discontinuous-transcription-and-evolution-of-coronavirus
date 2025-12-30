# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:47:47 2022

@author: DELL
"""
import numpy as np
from scipy import sparse
from collections import Counter
import re

def select_strc(ct_f, s, e, types = 'varna'):
    strc_selet = ct_f[s-1:e].copy()
    mask = (strc_selet['5'] < s) | (strc_selet['5'] > e)
    strc_selet['5'][mask] = 0
    if types == 'rna_strc':
        strc_selet['1'] = strc_selet['1'] - s +1
        strc_selet['3'] = strc_selet['3'] - s +1
        strc_selet['4'] = strc_selet['4'] - s +1
        strc_selet['5'][strc_selet['5']>0] = strc_selet['5'][strc_selet['5']>0] - s +1
        strc_selet['6'] = strc_selet['6'] - s +1
    return strc_selet

def ct2bedpe(dirs , ct_f, out_f, lens, types = 'ct'):    
    if types == 'ct':
        rna_strc = np.loadtxt(dirs + ct_f, usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
    if types == 'rna':
        rna_strc = np.loadtxt(dirs + ct_f, usecols=[0,1,2], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])
        
    outs = np.zeros(rna_strc.shape[0], dtype=[('c1','<U4'),('s1','<i4'),('e1','<i4'),('c2','<U4'),('s2','<i4'),('e2','<i4')])
    outs['c1'] = 'HS1'
    outs['c2'] = 'HS1'
    if types == 'ct':
        outs['s1'] = rna_strc['1']
        outs['s2'] = rna_strc['5']
        outs['e1'] = rna_strc['1'] + lens
        outs['e2'] = rna_strc['5'] + lens
    if types == 'rna':
        outs['s1'] = rna_strc['S']
        outs['s2'] = rna_strc['E']
        outs['e1'] = rna_strc['S'] + lens
        outs['e2'] = rna_strc['E'] + lens
        
    np.savetxt(dirs + out_f, outs[outs['s2']>0], fmt = '%s', delimiter='\t')

def ins2bedpe(dirs ,ins, out_f):    
    outs = np.zeros(ins.shape[0], dtype=[('c1','<U4'),('s1','<i4'),('e1','<i4'),('c2','<U4'),('s2','<i4'),('e2','<i4')])
    outs['c1'] = 'HS1'
    outs['c2'] = 'HS1'
    outs['s1'] = ins['1']
    outs['s2'] = ins['1'] 
    outs['e1'] = ins['1'] + ins['2']
    outs['e2'] = ins['1'] + ins['2']
    np.savetxt(dirs + out_f, outs, fmt = '%s', delimiter='\t', encoding='utf-8')

def del2bedpe(dirs , del_list, out_f): 
    del_list = np.array(del_list)
    outs = np.zeros(del_list.shape[0], dtype=[('c1','<U4'),('s1','<i4'),('e1','<i4'),('c2','<U4'),('s2','<i4'),('e2','<i4')])
    outs['c1'] = 'HS1'
    outs['c2'] = 'HS1'
    outs['s1'] = del_list[:,0]
    outs['s2'] = del_list[:,0]
    outs['e1'] = del_list[:,1]
    outs['e2'] = del_list[:,1]
    np.savetxt(dirs + out_f, outs, fmt = '%s', delimiter='\t', encoding='utf-8')



def Correlation_length(part):
    outs = []
    bound = [0]
    start = 0
    buff = part[0]
    for i,x in enumerate(part):
        if x != buff:
            outs.append( i - start )
            bound.append(i)
            buff = x
            start = i
    outs.append(i - start + 1)
    return np.array(outs),np.array(bound)
    
def get_matrix_npz(M1,M2):
    '''
    transform a three tuple into a matrix
    '''
    S1 = sparse.coo_matrix((M1['value'],(np.asarray(M1['S'], dtype=int),np.asarray(M1['E'], dtype= 'i4')))).toarray()
    S2 = sparse.coo_matrix((M2['value'],(np.asarray(M2['E'], dtype=int),np.asarray(M2['S'], dtype= 'i4')))).toarray()
    matrix = np.zeros((max(S2.shape),max(S1.shape)))
    matrix[np.nonzero(S1)] = S1[np.nonzero(S1)]
    matrix[np.nonzero(S2)] = S2[np.nonzero(S2)]
    return matrix

def get_oe(matrix):
    matrix = np.array(matrix).copy()
    lens = matrix.shape[0]
    mask = np.ones(lens)
    cut_off = lens * 0.0001
    matrix_num = matrix.copy()
    matrix_num[matrix_num > 0] = 1
    col_sum = matrix_num.sum(axis = 0)
    sites = np.arange(lens)
    for i in sites:
        if col_sum[i] <= cut_off:
            mask[i] = 0
    mask = mask == 1
    x_o = sites.copy()
    y_o = sites.copy()
    dis = np.zeros(lens)
    mask_o = mask.copy()
    exp_matrix = np.ones((lens, lens))
    x_e = sites.copy()
    y_e = sites.copy()
    for i in range(lens):
    ##obs_matrix
        diag_n = matrix[(x_o, y_o)]
        diag_line = diag_n[mask_o]
        if diag_line.shape[0] == 0:
            dis[i] = 0
        else :
            dis[i] = diag_line.mean()
            x_o = np.delete(x_o,-1)
            y_o = np.delete(y_o,0)
            mask_o = np.delete(mask_o,-1)
        
        ##exp_matrix
        if dis[i] == 0:
            x_e = np.delete(x_e,-1)
            y_e = np.delete(y_e,0)
            continue
        exp_matrix[(x_e, y_e)] = dis[i]
        exp_matrix[(y_e, x_e)] = dis[i]
        x_e = np.delete(x_e,-1)
        y_e = np.delete(y_e,0)

    oe = np.matrix(matrix/exp_matrix)
    return np.asarray(oe), dis

def apa(matrix, loop, s, e, lens, res, types = 'loop'):
    loop = np.array(loop, dtype = [('S','<i4'),('E','<i4'),('value','<f8')])
    ma_lens = lens*2 + 1
    ma = np.zeros([ma_lens, ma_lens])
    count = 0
    out_data = []
    for i in range(loop.shape[0]):
        start = int(loop['S'][i] / res )
        end = int(loop['E'][i] / res)
        if start <= lens or end >= (matrix.shape[0] - lens):
             continue
        if  (end - start) > s/res  and (end - start)  < e/res:
            apa_ma = matrix[(start - lens) :(start + lens + 1),(end - lens):(end + lens + 1)]
            # means = np.hstack([apa_ma[:lens,:lens],apa_ma[-lens:,-lens:]]).mean()
            
            ma += apa_ma
           # out_data.extend(list(matrix[(start-2):(start+3),(end-2):(end+3)].ravel()))
            if types == 'loop':
                out_data.append(apa_ma[lens-1:lens+2, lens-1:lens+2].mean())
            elif types == 'stem':
                out_data.append(apa_ma[np.arange(lens-lens,lens+lens+1)[1:-1],np.arange(lens-lens,lens+lens+1)[1:-1][::-1]].mean())
                out_data.append(apa_ma[np.arange(lens-lens+1,lens+lens+2)[1:-1],np.arange(lens-lens+1,lens+lens+2)[1:-1][::-1]].mean())
                out_data.append(apa_ma[np.arange(lens-lens-1,lens+lens)[1:-1],np.arange(lens-lens-1,lens+lens)[1:-1][::-1]].mean())
            elif types == 'diag':
                out_data.append( apa_ma[ np.arange(lens-lens,lens+lens+1)[1:-1], np.arange(lens-lens,lens+lens+1)[1:-1] ].mean() )
                out_data.append( apa_ma[ np.arange(lens-lens+1,lens+lens+2)[1:-1], np.arange(lens-lens+1,lens+lens+2)[1:-1] ].mean() )
                out_data.append( apa_ma[ np.arange(lens-lens-1,lens+lens)[1:-1], np.arange(lens-lens-1,lens+lens)[1:-1] ].mean() )
            count += 1
    return ma, count, out_data


def save_matrix(Ma,file_name):
    files = open(file_name,'w')
    for x in Ma:
        strs = 'MK\t'+str(int(x[0])*10) +'\t'+ str(int(x[0]+1)*10)+'\t'+'MK\t'+str(int(x[1])*10) +'\t'+ str(int(x[1]+1)*10)+'\t'+str(x[2])+'\n'
        files.write(strs)
    files.close()
    
def switch_ma(ma):
    ma = ma.copy()
    list_s = ma['S'].copy()
    list_e = ma['E'].copy()
    mask1 = ma['S'] > ma['E']
    ma['S'][mask1] = list_e[mask1]
    ma['E'][mask1] = list_s[mask1]
    return ma    

def get_com(matrix):
    matrix = matrix.copy()
    col_sum = matrix.sum(axis=0)
    matrix_correct = np.matrix(col_sum).T * np.matrix(col_sum) / col_sum.sum() #ma_sum#
    matrix = np.asarray(matrix) / np.asarray(matrix_correct)
    matrix[np.isnan(matrix)] = 0
    return matrix    

def extract_nanopore(file_dir, file_name, res):
    files_in = open(file_dir + file_name )
    junc_sites = []
    for line in files_in:
        lines = line.split('\t')  
        lens = len(lines)
        iters = int(lens/2)
        if lens == 2:
            continue
        
        if lens % 2 != 0:
            print('die out')
            break
        
        
        for i in range(iters-1):
            start = int(lines[i*2+1])
            end = int(lines[i*2+2])
            junc_sites.append((int(start/res), int(end/res)))
    junc_sites = Counter(junc_sites)
    ma_junc = np.asarray([(key[0],key[1], junc_sites[key]) for key in  junc_sites], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])
    return ma_junc

def read_fasta(file_dir, file_name):
    files_in = open(file_dir + file_name )
    dict_seq = {}
    for line in files_in:
        lines  = line.split('\n')[0]
        if lines == '':
            continue
        
        if re.match('^>', lines):
            key = lines
            
        else:
            dict_seq[key] = lines
    return dict_seq
    
def read_junction(file_dir, filename, chrs):
    junc = np.loadtxt(file_dir + filename, usecols=[0,1,2,3,4,5,6], dtype=[('chr_s','<U12'),('site_s','<i4'),('sign_s','<U1'),('chr_e','<U12'),('site_e','<i4'),('sign_e','<U1'),('direction','<i4')])
    junc = junc[junc['chr_s']!=junc['chr_e']]
    junc = junc[(junc['chr_s'] == chrs) | (junc['chr_e'] == chrs)]

    mask = junc['chr_e'] == chrs
    junc['chr_e'][mask] = junc['chr_s'][mask]
    junc['chr_s'][mask] = chrs
    tem_site = junc['site_s'][mask]
    junc['site_s'][mask] = junc['site_e'][mask]
    junc['site_e'][mask] = tem_site
    junc['direction'] = 1
    junc['direction'][mask] = -1
    return junc
    
def junction_to_bed(file_dir, filename, chrs):
    junc = np.loadtxt(file_dir + filename, usecols=[0,1,7,3,4], dtype=[('chr_s','<U12'),('site_s','<i4'),('site_s2','<i4'),('chr_e','<U12'),('site_e','<i4')])
    junc = junc[junc['chr_s']!=junc['chr_e']]
    junc = junc[(junc['chr_s'] == chrs) | (junc['chr_e'] == chrs)]

    mask = junc['chr_s'] == chrs
    junc['chr_s'][mask] = junc['chr_e'][mask]
    junc['chr_e'][mask] = chrs
    tem_site = junc['site_s'][mask]
    junc['site_s'][mask] = junc['site_e'][mask]
    junc['site_e'][mask] = tem_site
    return junc    

def split_ct(ct_in):
    st = np.zeros_like(ct_in['5'])
    st[ct_in['5'] > ct_in['1']] = 1
    st[(ct_in['5'] < ct_in['1']) & (ct_in['5'] !=0)] = -1
    st_type = []
    dict_st = {}
    sum_score = 0
    count = 1
    tem_score = 0 
    gap_s = 0
    bound = [0]
    for i,x in enumerate(st):
        sum_score += x
        if tem_score == 0 and sum_score == 0:
             st_type.append( 'gap_' + str(count) )
            
        elif tem_score == 0 and sum_score != 0:
             st_s = i
             st_type.append( 'stem-loop_' + str(count) )
             gap_e = i - 1
             dict_st['gap_' + str(count)] = (gap_s, gap_e)
            
        elif tem_score != 0 and sum_score != 0:
             st_type.append( 'stem-loop_' + str(count) )
            
        elif tem_score != 0 and sum_score == 0:    
             st_type.append( 'stem-loop_' + str(count) )
             st_e = i
             gap_s = i + 1
             dict_st['stem-loop_' + str(count)] = (st_s,st_e)  
             count += 1
             dict_st['gap_' + str(count)] = (gap_s, len(st) - 1) 
             bound.append(i)
        
            # if tem_score == 0:
            #     st_s = i
            #     gap_e = i - 1
            #     dict_st['gap_' + str(count)] = (gap_s, gap_e)    
        # elif sum_score > 0 :
        #    st_type.append( 'stem-loop_' + str(count) )
           
        tem_score = sum_score
    bound.append(i)
    return st_type, dict_st, bound

def extract_matrix(file_dir, file_name, res):
    data = np.loadtxt(file_dir + file_name, usecols=[0,1,2], dtype=[('S','<i4'),('E','<i4'),('value','<f8')])
    data['S'] = data['S'] / res
    data['E'] = data['E'] / res
    data['value'][data['value'] < 0] = 0
    data['value'][np.isnan(data['value'])] = 0
    mat = get_matrix_npz(data, data)
    oe, dis = get_oe(mat)
    return mat, oe, dis

def classify_rna_data(input_rna, dis, freq, binsize, discard):
    data = input_rna.copy()
    distance = data['E'] - data['S']
    rna_long = data[distance >= dis]
    rna_long_hf = rna_long[rna_long['value'] > freq]
    
    mask_c = (rna_long_hf['S'] <= discard[0]) & (rna_long_hf['value']>discard[1])
    
    hist_data = []
    s = 0
    e = binsize
    lens = data['E'].max()
    while True:
        distance_hf = rna_long_hf[~mask_c]['E'] - rna_long_hf[~mask_c]['S']
        mask = (distance_hf >= s) & (distance_hf <= e)
        hist_data.append(rna_long_hf[~mask_c]['value'][mask].sum())
        s = e
        e += binsize
        if e >= lens:
            break
    
    return rna_long_hf[mask_c], rna_long_hf[~mask_c], hist_data

def get_max_line(lines):
    lists = list(lines)
    out = []
    n = 0
    for l in lists:
        if l == '|':
            n+=1
        else:
            n = 0
        out.append(n)
    return max(out)

def get_shannon(matrix, n, res):
    (x,y) = np.diag_indices_from(matrix)
    shannon_ma = np.zeros_like(matrix)
    shan_out = []
    true_shan_out = []
    x = np.delete(x, 0)
    y = np.delete(y, -1)
    for i in np.arange(n-1):
        shannon_ma[x,y] = matrix[x,y].copy()
        shannon_ma[y,x] = matrix[y,x].copy()
        x = np.delete(x, 0)
        y = np.delete(y, -1)
    shannon_sum = shannon_ma.sum(axis=1)
    p_ma = shannon_ma / shannon_sum
    p_ma[np.isnan(p_ma)] = 0
    for i,arrs in enumerate(p_ma.T):
        p = arrs[arrs>0]
        shan = np.sum(p*np.log2(p)) * -1
        shan_out.append(shan)
        for iters in range(res):
            true_shan_out.append(shan)       
    return np.array(shan_out), np.array(true_shan_out)
