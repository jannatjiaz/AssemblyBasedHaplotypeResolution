import argparse
import pysam 
import pandas as pd
import random
import os 
import sys
from os import system
import re
import numpy as np
import multiprocessing as mp
import csv 
from intervaltree import Interval, IntervalTree
from math import floor
import random
from random import random 


#save reads
def save_list(MyList, output_file_name):
        MyFile=open(output_file_name,'w')
        MyList=map(lambda x:x+'\n', MyList)
        MyFile.writelines(MyList)
        MyFile.close()




parser = argparse.ArgumentParser(description='usage:  split_isoseq_reads.py  --sample sample --hap1_isoseq_reads hap1_isoseq_reads --hap2_isoseq_reads hap2_isoseq_reads --sample sample --hap1_contig_lengths hap1_contig_lengths --hap2_contig_lengths hap2_contig_lengths --haplotype1_read_dict haplotype1_read_dict --haplotype2_read_dict haplotype2_read_dict --outputdir outputdir')
parser.add_argument('--hap1_isoseq_reads', help='reads with hap1 snsp',required=True)
parser.add_argument('--hap2_isoseq_reads', help='reads with hap2 snps',required=True)
parser.add_argument('--haplotype1_read_dict', help='reads dictionary haplotype1',required=True)
parser.add_argument('--haplotype2_read_dict', help='reads dictionary haplotype2',required=True)
parser.add_argument('--sample', help='Input sample name',required=True)
parser.add_argument('--chromosome', help='Input chromosome',required=True)
parser.add_argument('--hap1_contig_lengths', help='contigs lengths file for hap1 allele',required=True)
parser.add_argument('--hap2_contig_lengths', help='contigs lengths file for hap2 allele',required=True)
parser.add_argument('--outputdir', help='output directory',required=True)
args = parser.parse_args()

chromosome=args.chromosome
sample=args.sample
hap1_isoseq_readsfile = args.hap1_isoseq_reads
hap2_isoseq_readsfile = args.hap2_isoseq_reads
hap1_contig_lengthsfile = args.hap1_contig_lengths
hap2_contig_lengthsfile= args.hap2_contig_lengths
haplotype1_read_dict_file = args.haplotype1_read_dict
haplotype2_read_dict_file = args.haplotype2_read_dict
outputdir=args.outputdir

print(sample)

hap1_read_dict = np.load("{}".format(haplotype1_read_dict_file)).item()
hap2_read_dict = np.load("{}".format(haplotype2_read_dict_file)).item()

#load in phased reads 
hap1_isoseq_reads = list(pd.read_csv("{}".format(hap1_isoseq_readsfile), sep=' ', header=None)[0])
hap2_isoseq_reads = list(pd.read_csv("{}".format(hap2_isoseq_readsfile), sep=' ', header=None)[0])
#turn phased reads into a dictionaty - its faster
hap1_isoseq_reads_dict={}
for read in hap1_isoseq_reads:
        hap1_isoseq_reads_dict[read]=1

#turn phased reads into a dictionaty - its faster
hap2_isoseq_reads_dict={}
for read in hap2_isoseq_reads:
        hap2_isoseq_reads_dict[read]=1

hap1_input  = list(hap1_read_dict.keys())
hap1_size = int(round(len(hap1_input)/19))
hap1_chunks = [hap1_input[i:i+int(hap1_size)] for i in range(0, len(hap1_input), int(hap1_size))]  

#make a list of wildtype contigs
hap1_contig_list = list(pd.read_csv(hap1_contig_lengthsfile, sep='\t', header=None)[0])
hap2_contig_list = list(pd.read_csv(hap2_contig_lengthsfile, sep='\t', header=None)[0])

hap2_contig_len = pd.read_csv(hap2_contig_lengthsfile, sep='\t', header=None)
hap2_contig_len_dict={}
for item in range(len(list(hap2_contig_len[0]))):
	hap2_contig_len_dict[hap2_contig_len[0][item]]=int(hap2_contig_len[1][item])


hap1_contig_len = pd.read_csv(hap1_contig_lengthsfile, sep='\t', header=None)
hap1_contig_len_dict={}
for item in range(len(list(hap1_contig_len[0]))):
	hap1_contig_len_dict[hap1_contig_len[0][item]]=int(hap1_contig_len[1][item])

poor_quality_reads = []
equal_mapping_reads=[]
equal_mapping_reads_pos=[]
hap2_reads=[]
hap1_reads=[]
def split_reads_parallel(hap1_chunks):
        chunk = hap1_chunks
        print("sent to core")
        poor_quality_reads = []
        equal_mapping_reads_pos=[]
        equal_mapping_reads=[]
        hap2_reads=[]
        hap1_reads=[]
        for hap1_read in chunk:
                if hap1_read in hap2_read_dict:
                #if the hetreozygous snps indicate that the reads should be in a specific group:
                        if hap1_read in hap1_isoseq_reads_dict and hap1_read not in hap2_isoseq_reads_dict:
                                hap1_reads.append(hap1_read)
                        elif hap1_read in hap2_isoseq_reads_dict and hap1_read not in hap1_isoseq_reads_dict:
                                hap2_reads.append(hap1_read)
                        #if all reads are mapped
                        elif hap1_read_dict[hap1_read]['mapping_read']=='mapped' and hap2_read_dict[hap1_read]['mapping_read']=='mapped':
                                #if all the mapping is bad 
                                if int(hap1_read_dict[hap1_read]['read_mapq'])<10 and int(hap2_read_dict[hap1_read]['read_mapq'])<10:
                                        poor_quality_reads.append(hap1_read)
                                #if the wildtpye mapping is bad but the hap1 is not 
                                elif int(hap1_read_dict[hap1_read]['read_mapq'])>10 and int(hap2_read_dict[hap1_read]['read_mapq'])<10:
                                        hap1_reads.append(hap1_read)
                                #if the hap1 mapping is bad but the hap2 is not 
                                elif int(hap1_read_dict[hap1_read]['read_mapq'])<10 and int(hap2_read_dict[hap1_read]['read_mapq'])>10:
                                        hap2_reads.append(hap1_read)
                                #if the hap1 read mapping is better than the wildtype read mapping and none of the hap1 mapping is 0
                                elif int(hap1_read_dict[hap1_read]['read_mapq']) > int(hap2_read_dict[hap1_read]['read_mapq']):
                                        hap1_reads.append(hap1_read)
                                #if the wildtype read mapping is better than the hap1 read mapping and none of the hap1 mapping is 0
                                elif int(hap1_read_dict[hap1_read]['read_mapq']) < int(hap2_read_dict[hap1_read]['read_mapq']):
                                        hap2_reads.append(hap1_read)
                                #if the mapping is the same 
                                elif int(hap1_read_dict[hap1_read]['read_mapq']) == int(hap2_read_dict[hap1_read]['read_mapq']):
                                        equal_mapping_reads.append(hap1_read)
                                        equal_mapping_reads_pos.append(hap1_read_dict[hap1_read]['read_pos'])
                                else:
                                        print('you missed something in mapq', hap1_read_dict[hap1_read]['read_mapq'], hap2_read_dict[hap1_read]['read_mapq'])
                        #if only the hap1 reads map - dont worry about mapq
                        elif (hap1_read_dict[hap1_read]['mapping_read']=='mapped') and (hap2_read_dict[hap1_read]['mapping_read']=='unmapped'):
                                hap1_reads.append(hap1_read)
                        #if only the hap2_reads map - dont worry about mapq
                        elif (hap2_read_dict[hap1_read]['mapping_read']=='mapped') and (hap1_read_dict[hap1_read]['mapping_read']=='unmapped'):
                                hap2_reads.append(hap1_read)
                        #no mapping for any of the read
                        elif hap1_read_dict[hap1_read]['mapping_read']=='unmapped' and hap2_read_dict[hap1_read]['mapping_read']=='unmapped':
                                poor_quality_reads.append(hap1_read)
                        else:
                                print("you missed something ", hap1_read )
                                poor_quality_reads.append(hap1_read)
                else:
                        #read missing in hap2 dict
                        hap1_reads.append(hap1_read)
        
        return poor_quality_reads, hap2_reads, hap1_reads, equal_mapping_reads, equal_mapping_reads_pos





pool = mp.Pool(9) #change if you change the number of cores the job is submitted to
output = pool.map(split_reads_parallel, hap1_chunks)
pool.close()
pool.join()

Array_Output = np.asarray(output)
poor_quality_reads = np.concatenate(Array_Output[:,0])
haplotype2_reads = np.concatenate(Array_Output[:,1])
haplotype1_reads = np.concatenate(Array_Output[:,2])
equal_mapping_reads = np.concatenate(Array_Output[:,3])
equal_mapping_reads_pos = np.concatenate(Array_Output[:,4])

haplotype2_reads=haplotype2_reads.tolist()
haplotype1_reads=haplotype1_reads.tolist()
equal_mapping_reads=equal_mapping_reads.tolist()
equal_mapping_reads_pos=equal_mapping_reads_pos.tolist()

equal_mapping_reads_df = pd.DataFrame({"read":equal_mapping_reads, "reads_pos":equal_mapping_reads_pos})
equal_mapping_reads_df = equal_mapping_reads_df.sort_values(by=['reads_pos'])
equal_mapping_reads_df = equal_mapping_reads_df.reset_index(drop=True)

for haplotype1_read in hap1_read_dict:
        if haplotype1_read not in hap2_read_dict:
                haplotype1_reads.append(haplotype1_read)

#save only the reads which are actually assigned 
save_list(haplotype1_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype1_reads_assigned_',chromosome,'_',sample,'.txt'))
save_list(haplotype2_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype2_reads_assigned_',chromosome,'_',sample,'.txt'))


#append the equal mapping reads to haplotype2 reads
for read in range(0,len(equal_mapping_reads_df),2):
        haplotype2_reads.append(equal_mapping_reads_df["read"][read])


#append the equal mapping reads to haplotype2 reads
for read in range(1,len(equal_mapping_reads_df),2):
        haplotype1_reads.append(equal_mapping_reads_df["read"][read])

save_list(poor_quality_reads, "{}{}{}{}{}{}".format(outputdir,'poor_quality_reads_',chromosome,'_',sample,'.txt'))
save_list(equal_mapping_reads, "{}{}{}{}{}{}".format(outputdir,'equal_mapping_reads_',chromosome,'_',sample,'.txt'))
save_list(haplotype1_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype1_reads_',chromosome,'_',sample,'.txt'))
save_list(haplotype2_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype2_reads_',chromosome,'_',sample,'.txt'))


