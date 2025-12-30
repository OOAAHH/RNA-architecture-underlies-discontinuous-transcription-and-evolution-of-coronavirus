# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 15:09:30 2023

@author: Administrator
"""
import numpy as np
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
color = ['lightgreen','lightgreen', 'pink', 'lightblue', 'pink', 'lightblue', 'pink', 'lightblue', 'pink', 'lightblue','pink', 'lightblue','pink', 'lightblue']

def plot_matrix_ct(matrix, pcs, ct, title, s , e, r):
    '''
    plot a compartment and heatmap mix fig
    '''
    ##
    mean = [x.mean() for x in pcs]
    pcs = [x[s:e] for x in pcs]

    ##ax size
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    size_heatmap = [left, bottom, width, height]
    size_colorbar = [left + width + 0.04, bottom, width/20, height]
    fig = plt.figure(figsize=(12,12))
    
    ##ct
    ct = ct[s:e]
    mask_ct = (ct['5'] < (e-1)) & (ct['5'] > s)
    ct = ct[mask_ct]
    start = s
    col_ric = (ct['1']- start - 1) / r
    raw_ric = (ct['5'] - start - 1) / r
    
    ##heatmap 
    s = int(s / r)
    e = int(e / r)
    matrix = matrix[s:e, s:e]
    ax = fig.add_axes(size_heatmap)
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero,95)
    sc = ax.imshow(matrix, cmap = plt.get_cmap('bwr'), vmax = vmax, vmin = 0-vmax,aspect = 'auto', 
                   interpolation = 'none', origin = 'lower')
    
    plt.scatter(raw_ric[col_ric < raw_ric],col_ric[col_ric < raw_ric], s = r,marker = 's', color = 'none', edgecolors ='black')
    plt.xticks(np.linspace(0, e - s, 5),np.linspace(s*r, e*r,5, dtype='i4'))
    plt.yticks(np.linspace(0, e - s, 5),np.linspace(s*r, e*r,5, dtype='i4'))
    ax.set_xlabel('Bp')
    ax.set_ylabel(title[-1],fontsize=20)
    ax = fig.add_axes(size_colorbar)
    fig.colorbar(sc,cax = ax)
    
    ##compartments
    lens = len(pcs)
    stept_h = 0.2 / lens
    bottom_h = bottom + height
    for i in range(lens):
        ax = fig.add_axes([left, bottom_h, width, stept_h])
        ax.fill_between(np.arange(len(pcs[i])), pcs[i])
        ax.set_xlim((0,len(pcs[i])-1))
        ax.set_xticks([])
        ax.set_yticks([mean[i]])
        ax.set_ylabel(title[i], fontsize=10, rotation = 'horizontal', labelpad = 20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        bottom_h = bottom_h + stept_h
    
    return fig

def plot_matrix_ct2(matrix, pcs, ct, title, s , e, r):
    '''
    plot a compartment and heatmap mix fig
    '''
    ##
    mean = [x.mean() for x in pcs]
    pcs = [x[s:e] for x in pcs]

    ##ax size
    left, width = 0.15, 0.65
    bottom, height = 0.1, 0.65
    size_heatmap = [left, bottom, width, height]
    size_colorbar = [left + width + 0.04, bottom, width/20, height]
    fig = plt.figure(figsize=(8,8))
    
    ##ct
    ct = ct[s:e]
    mask_ct = (ct['5'] < (e-1)) & (ct['5'] > s)
    ct = ct[mask_ct]
    start = s
    col_ric = (ct['1']- start - 1) / r
    raw_ric = (ct['5'] - start - 1) / r
    
    ##heatmap 
    s = int(s / r)
    e = int(e / r)
    matrix = matrix[s:e, s:e]
    ax = fig.add_axes(size_heatmap)
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero,95)
    sc = ax.imshow(matrix, cmap = plt.get_cmap('bwr'), vmax = vmax, vmin = 0-vmax,aspect = 'auto', 
                   interpolation = 'none', origin = 'lower')
    
    plt.scatter(raw_ric[col_ric < raw_ric],col_ric[col_ric < raw_ric], s = r,marker = 's', color = 'none', edgecolors ='black')
    plt.plot([50/r,50/r], [50/r, (e-s)-50/r], linestyle = ':', color = 'blue', alpha = 0.25)
    plt.plot([50/r, (e-s)-50/r], [(e-s)-50/r, (e-s)-50/r], linestyle = ':', color = 'blue', alpha = 0.25)
    plt.xticks(np.linspace(0, e - s, 5),np.linspace(s*r, e*r,5, dtype='i4'))
    plt.yticks(np.linspace(0, e - s, 5),np.linspace(s*r, e*r,5, dtype='i4'))
    ax.set_xlabel('Bp')
    ax.set_ylabel(title[-1],fontsize=20)
    ax = fig.add_axes(size_colorbar)
    fig.colorbar(sc,cax = ax)
    
    ##compartments
    lens = len(pcs)
    stept_h = 0.2 / lens
    bottom_h = bottom + height
    for i in range(lens):
        ax = fig.add_axes([left, bottom_h, width, stept_h])
        ax.fill_between(np.arange(len(pcs[i])), pcs[i])
        ax.set_xlim((0,len(pcs[i])-1))
        ax.set_xticks([])
        ax.set_yticks([mean[i]])
        ax.set_ylabel(title[i], fontsize=10, rotation = 'horizontal', labelpad = 20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        bottom_h = bottom_h + stept_h
    
    return fig

def joint_heatmap_genome(matrix, dict_gene, title, r, vmax, cmap):
    '''
    plot a compartment and heatmap mix fig
    '''
    ##ax size
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    size_heatmap = [left, bottom, width, height]
    size_colorbar = [0.8, 0.775, 0.015, 0.15]
    n = int(50 / r)
    fig = plt.figure(figsize=(10,10))

    ##heatmap 
    ax = fig.add_axes(size_heatmap)  
    shp = matrix.shape
    lenth = shp[0] * r 
    ma = np.log2(matrix+1)
    sc = ax.imshow(ma, cmap = plt.get_cmap(cmap), vmax = vmax, vmin = 0
              , aspect = 'auto',extent = (0, lenth, lenth, 0))
    ax.spines[['top','right']].set_visible(False)    
    ax.set_xlabel('Bp')
    ax.set_ylabel(title,fontsize=20)
    ax = fig.add_axes(size_colorbar)
    fig.colorbar(sc,cax = ax)
      
    ## genome
    ax = fig.add_axes([left, bottom + height + 0.01, width, 0.05])
    i = 0
    for key in dict_gene:
        s,e = dict_gene[key]
        s = s / r
        e = e / r
            
        if key == 'leader':
            plt.fill_between([s,e],[0.25,0.25], -0.25, color = 'black', alpha = 0.5)
       
        elif key == 'UTR_5':
            plt.plot([s,e],[0,0], lw = 2, color = 'grey')
        
        elif key == 'UTR_3':
            plt.fill_between([s,e],[0.25,0.25], -0.25, color = 'grey', alpha = 0.5)
        
        else :
            if i % 2 == 0:
                plt.fill_between([s,e],[0.2,0.2], -0.1, color = color[i], alpha = 0.5)
            if i % 2 == 1:
                plt.fill_between([s,e],[0.1,0.1], -0.2, color = color[i], alpha = 0.5)
            i += 1
            
    plt.ylim(-0.25,0.3)
    plt.yticks([])
    genome_lens = int(dict_gene['UTR_3'][1] / r) + 1
    plt.xlim(0, genome_lens)

    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.spines[['left', 'right', 'bottom']].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    ax.tick_params(which='major', width=1.00, length=5)
    plt.xticks(np.arange(0, genome_lens, 100 * n), np.arange(0, genome_lens, 100 * n) * r)
    
    ## genome

    ax = fig.add_axes([left + width + 0.01, bottom, 0.05,  height])
    i = 0
    for key in dict_gene:
        s,e = dict_gene[key]
        s = s / r
        e = e / r
        s = genome_lens - s
        e = genome_lens - e
        if key == 'leader':
            plt.fill_betweenx([s,e],[0.25,0.25], -0.25, color = 'black', alpha = 0.5)
       
        elif key == 'UTR_5':
            plt.plot([0,0], [s,e], lw = 2, color = 'grey')
        
        elif key == 'UTR_3':
            plt.fill_betweenx([s,e],[0.25,0.25], -0.25, color = 'grey', alpha = 0.5)
        
        else :
            if i % 2 == 0:
                plt.fill_betweenx([s,e],[0.2,0.2], -0.1, color = color[i], alpha = 0.5)
            if i % 2 == 1:
                plt.fill_betweenx([s,e],[0.1,0.1], -0.2, color = color[i], alpha = 0.5)
            i += 1
            
    plt.xlim(-0.25,0.3)
    plt.yticks([])
    plt.xticks([])
    plt.ylim(0, genome_lens)

    ax.tick_params(right=True, labelright=True, left=False, labelleft=False)
    ax.spines[['left', 'top', 'bottom']].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    ax.tick_params(which='major', width=1.00, length=5)
    plt.yticks(np.arange(genome_lens%(100*n), genome_lens+100*n, 100*n), np.arange(0, genome_lens, 100*n)[::-1] * r)
    
    return fig

def plot_tri_matrix(matrix, ct_in, title, S, E, H, r, interp):
    '''
    plot a compartment and heatmap mix fig
    '''
    ##ax size
    matrix[np.isnan(matrix)] = 0
    matrix = np.log2(matrix + 1)
    ratio = (H / (E - S))
    left, width = 0.1, 0.65
    if ratio > 0.35 :
        ratio = 0.35
    bottom, height = 1 -  0.65 * ratio  - 0.1, 0.65 * ratio
    size_heatmap = [left, bottom, width, height]
    size_colorbar = [left + width + 0.04, bottom, width/160, height/2]
    fig = plt.figure(figsize=(16 * ratio * 2**0.5 , 12))
    
    ##ct 
    s = S * r
    e = E * r
    ct = ct_in[s:e]
    mask_ct = (ct['5'] < (e-1)) & (ct['5'] > s)
    ct = ct[mask_ct]
    start = s
    col_ric = (ct['1']- start - 1) / r
    raw_ric = (ct['5'] - start - 1) / r
    
    ##heatmap 
    ax = fig.add_axes(size_heatmap)
    shp = matrix.shape
    tri_ma = np.zeros([H, E-S])
    x,y = np.diag_indices_from(matrix)
    for i in np.arange(H):
        if (S - i) > 0 and (E + i) < shp[0]:
            tri_ma[H-i-1] = matrix[x[S+i:E+i],y[S-i:E-i]]
            
        elif (S - i) < 0:
            step = i - S
            tri_ma[H-i-1][step:] = matrix[x[S+i+step:E+i],y[S-i+step:E-i]]
            
        elif (E + i) > shp[0]:
            step = shp[0] - (E+i)
            tri_ma[H-i-1][:step] = matrix[x[S+i:E+i+ step],y[S-i:E-i + step]]
        
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero,85)
    
    sc = ax.imshow(tri_ma, vmax = vmax, vmin = 0, aspect = 'auto', interpolation = interp, cmap = plt.get_cmap('pink_r'))
    ax.set_ylabel(title[-1],fontsize=20)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    plt.xlim([0,E-S])
    plt.xticks(np.arange(0, E-S + 1, r * 10), np.arange(S, E+1, r*10) * r)
    plt.yticks([])
    ylim = plt.ylim()
    ax.spines[['right', 'left', 'bottom']].set_visible(False)
    
    ax = fig.add_axes(size_colorbar)
    fig.colorbar(sc,cax = ax)
    
    ##scatter stem loop  
    stept_h = height
    bottom_h = bottom - stept_h 
    ax = fig.add_axes([left, bottom_h, width, stept_h])
    ax.imshow(tri_ma[::-1], vmax = vmax, vmin = 0, aspect = 'auto', interpolation = interp, cmap = plt.get_cmap('pink_r'))
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
    plt.scatter((raw_ric+col_ric)/2, (col_ric-raw_ric) /2, s = r, marker = 's', color = 'none', edgecolors ='blue', alpha=0.05)
    plt.xlim([0,E-S])
    plt.ylim(ylim)
    plt.yticks([])
    plt.xticks([])
    
    ##line stem loop 
    bottom_h = bottom_h - stept_h
    ax = fig.add_axes([left, bottom_h, width, stept_h])
    for seq in zip(raw_ric, col_ric):
        s = seq[0] 
        e = seq[1] 
        radius = (e - s) / 2
        a,b = ((s+e)/2, 0)
        x = np.arange(a - radius, a + radius + 1, 0.1 )
        y = b + np.sqrt(radius ** 2 - (x - a) ** 2)
        plt.plot(x , y, lw = 0.15)  
        plt.plot(x , -y, lw = 0.15)  
    ax.spines[['top', 'right', 'left']].set_visible(False)
    plt.xlim([0,E-S])
    plt.ylim([0,H])
    ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    plt.xticks(np.arange(0, E-S + 1, r * 10), np.arange(S, E+1, r*10) * r)
    plt.yticks([])
    ax.set_xlabel('Bp')
    return fig

def plot_tri_matrix2(matrix, ct_in, title, s, e, r, interp, cmap):
    '''
    plot a compartment and heatmap mix fig
    '''
    ##ax size
    S = int(s / r)
    E = int(e / r)
    matrix = matrix.copy()
    matrix[np.isnan(matrix)] = 0
    matrix = np.log2(matrix + 1)
    H = int((E - S)/2)
    left, width = 0.1, 0.75
    ratio = 0.5
    bottom, height = 1 -  0.75 * ratio , 0.75 * ratio
    size_stem = [left, bottom - 0.1, width, height]
    stept_h = height
    bottom_h = bottom - stept_h - 0.04
    size_heatmap = [left, bottom_h - 0.1, width, stept_h]
    size_colorbar = [left + width - 0.1, bottom_h- 0.05, width/40, height/2]
    fig = plt.figure(figsize=(12 * ratio * 2**0.5 , 12))
    
    ##ct 
    ct = ct_in[s:e]
    mask_ct = (ct['5'] < (e-1)) & (ct['5'] > s)
    ct = ct[mask_ct]
    start = s
    col_ric = (ct['1']- start - 1) / r
    raw_ric = (ct['5'] - start - 1) / r
    
    ##heatmap 
    ax = fig.add_axes(size_heatmap)
    shp = matrix.shape
    ma = np.zeros_like(matrix) * np.nan
    ma[S:E,S:E] = matrix[S:E,S:E]
    tri_ma = np.zeros([H, E-S]) * np.nan
    x,y = np.diag_indices_from(matrix)
    for i in np.arange(H):
        if (S - i) >= 0 and (E + i) <= shp[0]:
            tri_ma[H-i-1] = ma[x[S+i:E+i],y[S-i:E-i]]
                  
        elif (S - i) < 0:
            step = i - S
            tri_ma[H-i-1][step:] = ma[x[S+i+step:E+i],y[S-i+step:E-i]] 
            
        elif (E + i) > shp[0]:
            step = shp[0] - (E+i)
            tri_ma[H-i-1][:step] = ma[x[S+i:E+i+ step],y[S-i:E-i + step]] 
        
    nonzero = matrix[S:E,S:E][np.nonzero(matrix[S:E,S:E])]
    vmax = np.percentile(nonzero,85)
    
    ##scatter stem loop 
    sc = ax.imshow(tri_ma[::-1], vmax = vmax, vmin = 0, aspect = 'auto'
                   , interpolation = interp, cmap = plt.get_cmap(cmap))
    plt.scatter((raw_ric+col_ric)/2, (col_ric-raw_ric) /2 , s = 50, marker = 'v',
                color = 'none', edgecolors ='red', alpha=0.25)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    plt.ylim(H,0)
    plt.xlim([0,E-S])
    plt.xticks(np.arange(0, E-S + 1, r * 10), np.arange(S, E+1, r*10) * r)
    plt.yticks([])
    ax.spines[['right', 'left', 'bottom']].set_visible(False)
    
    ax = fig.add_axes(size_colorbar)
    fig.colorbar(sc,cax = ax)
    
    ##line stem loop  

    ax = fig.add_axes(size_stem)
   
    for seq in zip(raw_ric, col_ric):
        s = seq[0] 
        e = seq[1] 
        radius = (e - s) / 2
        a,b = ((s+e)/2, 0)
        x = np.arange(a - radius, a + radius + 1, 0.01 )
        y = b + np.sqrt(radius ** 2 - (x - a) ** 2)
        plt.plot(x , y, lw = 0.15, color = 'grey')  
        plt.plot(x , -y, lw = 0.15, color = 'grey')  
    ax.spines[['top', 'right', 'left']].set_visible(False)
    plt.xlim([0,E-S])
    plt.ylim([0,H])
    ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    #plt.xticks(np.arange(0, E-S + 1, r * 10), np.arange(S, E+1, r*10) * r)
    plt.xticks(np.arange(0, E-S + 1, r * 10), (' '*np.arange(S, E+1, r*10).shape[0]).split(' ')[:-1])
    # plt.xticks([])
    plt.yticks([])
    ax.set_title(title,fontsize=20)
    #ax.set_xlabel('Bp')
    return 

