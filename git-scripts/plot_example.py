# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 20:31:26 2023

@author: Administrator
"""

from Coronavirus_ct_shannon_data import strc_ric_pedv_cell, strc_ric_pedv_virion, strc_ric_sars2_virion, strc_ric_pdcov_virion
from Coronavirus_ct_shannon_data import shannon_ric_pedv_virion, shannon_ric_sars2_virion, shannon_ric_pdcov_virion
from Coronavirus_hic_data import ma_mers_com_cell_s, ma_sars_com_cell, ma_sars_ric_virion, ma_sars_rna_cell, ma_pedv_ric_cell_kr, ma_pedv_ric_virion_kr, ma_sars_splash_cell, ma_ric_pdcov_virion
from plots import plot_tri_matrix2, plot_matrix_ct, plot_matrix_ct2
from PEDV_deletion import pedv_long_del
from SARS2_deletion import sars2_long_del
from PDCoV_deletion import pdcov_long_del
from MERS_deletion import mers_long_del
import seaborn as sns
from genome_382_short import dict_gene_recom
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sns.set_theme(style="ticks")
color = ['lightgreen','lightgreen', 'pink', 'lightblue', 'pink', 'pink', 'lightblue', 'lightblue', 'lightblue', 'lightblue', 'pink', 'lightblue', 'lightblue', 'lightblue']
cmap_list = ['afmhot','viridis', 'cividis', 'RdYlGn', 'pink_r']
interpolation_list = ['antialiased', 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos', 'blackman']

plot_tri_matrix2(ma_sars_com_cell, strc_ric_sars2_virion, 'SARS-CoV-2:27250-29500Bp', 27250, 29500, 5, 'nearest', 'viridis')
plot_tri_matrix2(ma_pedv_ric_cell_kr, strc_ric_pedv_cell, 'PEDV:25000-25500Bp', 25000, 25500, 5, 'nearest', 'viridis')

# plot_matrix_ct(ma_sars_ric_virion, [shannon_ric_sars2_virion],strc_ric_sars2_virion, 
#                ['Shannon','SARS-CoV2:27250-29500Bp'], 27250, 29500, 5)
# plot_matrix_ct(ma_sars_com_cell, [shannon_ric_sars2_virion],strc_ric_sars2_virion, 
#                ['Shannon','SARS-CoV2:27250-29500Bp'], 27250, 29500, 5)
# plot_matrix_ct(ma_sars_rna_cell, [shannon_ric_sars2_virion],strc_ric_sars2_virion, 
#                ['Shannon','SARS-CoV2:27250-29500Bp'], 27250, 29500, 5)

# plot_tri_matrix2(ma_pedv_ric_virion, strc_ric_pedv_virion, 'PEDV deletion 20723-21322', 20650, 21400, 5, 'nearest', 'viridis')

# pp = PdfPages('PEDV-deletion-ric-cell-short-20bp.pdf')
# for pedv_del in sorted(pedv_long_del):
#     s = pedv_del[0]
#     e = pedv_del[1]
#     pp.savefig(plot_matrix_ct2(ma_pedv_ric_cell_kr, [shannon_ric_pedv_virion],strc_ric_pedv_virion, 
#                 ['Shannon','PEDV deletion:'+str(s) + '-' +str(e)], s-50, e+50, 5))
#     plt.close()
# pp.close()   

# pp = PdfPages('SARS-CoV-2-deletion-com.pdf')
# for sars_del in sorted(sars2_long_del)[:]:
#     s = sars_del[0]
#     e = sars_del[1]
#     pp.savefig(plot_matrix_ct2(ma_sars_com_cell, [shannon_ric_sars2_virion],strc_ric_sars2_virion, 
#                 ['Shannon','SARS-CoV2 deletion:'+str(s) + '-' +str(e)], s-50, e+50, 5))
#     plt.close()
# pp.close()

# pp = PdfPages('PDCoV-deletion-ric.pdf')
# for pd_del in sorted(pdcov_long_del)[:]:
#     s = pd_del[0]
#     e = pd_del[1]
#     pp.savefig(plot_matrix_ct2(ma_ric_pdcov_virion, [shannon_ric_pdcov_virion], strc_ric_pdcov_virion, 
#                 ['Shannon','PDCoV deletion:'+str(s) + '-' +str(e)], s-50, e+50, 5))
#     plt.close()
# pp.close()

# pp = PdfPages('MERS-deletion-ric.pdf')
# for ms_del in sorted(mers_long_del)[:]:
#     s = ms_del[0]
#     e = ms_del[1]
#     pp.savefig(plot_matrix_ct2(ma_mers_com_cell_s, [shannon_ric_pdcov_virion], strc_ric_pdcov_virion, 
#                 ['Shannon','MERS deletion:'+str(s) + '-' +str(e)], s-50, e+50, 5))
#     plt.close()
# pp.close()

# fig = plt.figure(figsize=(14,2))
# res = 1000
# i = 7
# for key in dict_gene_recom:
#     s,e = dict_gene_recom[key]
#     s = s / res
#     e = e / res
#     if key.split('_')[0] == 'gap':
#         plt.plot([s,e],[0,0], color = 'black', alpha = 0.5)
    
#     else:
#         plt.fill_between([s,e],[0.05,0.05], -0.05, color = color[i], alpha = 0.5)
#         i += 1
#     plt.xlim(27.5, 29.6)
#     plt.ylim(-0.25, 0.35)
# plt.yticks([])
# ax = plt.gca()
# ax.spines[['left', 'right', 'top']].set_visible(False)