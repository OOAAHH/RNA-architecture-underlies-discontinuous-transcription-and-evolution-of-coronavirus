# -*- coding: utf-8 -*-
# """
# Created on 2023.8.7
# Nanopore data analysis protocol
# @author: Zi

import os,argparse

def Nanoporepipeline(refseq, fqdata, v_chr, title):
    os.system("mkdir %s"%(title))
    os.system("minimap2 -d %s.min %s"%(refseq,refseq))#索引
    os.system("minimap2 -ax splice %s %s | samtools view -hF4 -b | samtools sort - > ./%s/%s_mapped.bam"%(refseq, fqdata, title, title))#比对以及过滤
    os.system("samtools view -h -o ./%s/%s_mapped.sam ./%s/%s_mapped.bam "%(title, title, title, title))#bam to sam
    os.system("python /store/zwen/workspace/Nanopore/npTranscript/sam_to_bed_linux.py %s ./%s/%s_mapped.sam ./%s/%s_mapped.bed "%(v_chr, title, title, title, title))#提取CIGAR信息
    os.system("python /store/zwen/workspace/Nanopore/npTranscript/counts_gene_linux.py %s ./%s/%s_mapped.bed ./%s/%s"%(v_chr, title, title, title, title))#计算基因表达


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="vRICseq data analysis pipeline-RNAseq. First-writtern by Max.")
    parser.add_argument("-r", "--refseq",type=str,required=True,help="Please type in the file of refseq.")
    parser.add_argument("-f","--fqdata",type=str,required=True,help="Please type in the file of fq sequencing data")
    parser.add_argument("-v","--v_chr",type=str,required=True,help="Please type in the chr of virus")
    parser.add_argument("-t", "--title", type=str, required=True,help="Please type in the name of outputfile;for example:SARSCOV2-1")
    Args = parser.parse_args()
    Nanoporepipeline(os.path.abspath(Args.refseq), os.path.abspath(Args.fqdata), Args.v_chr, Args.title)
    
    