# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:12:36 2022

@author: DELL
"""

import numpy as np
import re
from sys import argv
pyname, chrs, file_in_name, file_out_name = argv 

files_in = open(file_in_name)
file_out = open(file_out_name,'w')


for line in files_in:
     lines = line.split('\t')
     if len(lines) <= 4:
         continue
     
     if lines[2] == chrs:
         
         start = int(lines[3])
         # if start > 125 :
         #     continue
         info = lines[5]
         number = np.asarray(re.findall('\d+', info), dtype='i4')
         strs = np.asarray(re.findall('[A-z]', info))
         bed = []
         end = start
         
         for i,s in enumerate(strs):
             lens = number[i]
             if s == 'M' or s == 'D' :
                 end +=  lens
           
             elif s == 'N' :
                if lens <= 10:
                    end += lens
                    
                else:
                    bed.extend([str(start), str(end)])
                    start = end + lens
                    end = start
                    
         bed.extend([str(start), str(end)])
         file_out.writelines('\t'.join(bed)+'\n')
file_out.close()
files_in.close()




    





