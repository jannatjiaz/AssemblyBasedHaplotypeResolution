
import argparse
import pysam 
import pandas as pd
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
import argparse 


parser = argparse.ArgumentParser(description='usage:  compare_alignments.py  \
    --chromosome chromosome \
    --haplotype1_read_dict haplotype1_read_dict\
    --haplotype2_read_dict haplotype2_read_dict\
    --haplotype1_contig_lenghts haplotype1_contig_lenghts \
    --haplotype2_contig_lenghts haplotype2_contig_lenghts \
    --haplotype1_hic_snp_phased_reads haplotype1_hic_snp_phased_reads \
    --haplotype2_hic_snp_phased_reads haplotype2_hic_snp_phased_reads \
    --outputdir path')
parser.add_argument('--chromosome', help='chromosome name',required=True)
parser.add_argument('--haplotype1_read_dict', help='output of create_haplotype1_dict.py',required=True)
parser.add_argument('--haplotype2_read_dict', help='output of create_haplotype2_dict.py',required=True)
parser.add_argument('--haplotype1_hic_snp_phased_reads', help='hic reads that contain informative snps for haplotype1',required=True)
parser.add_argument('--haplotype2_hic_snp_phased_reads', help='hic reads that contain informative snps for haplotype 2',required=True)
parser.add_argument('--haplotype1_contig_lenghts', help='file of contig lengths',required=True)
parser.add_argument('--haplotype2_contig_lenghts', help='file of contig lengths',required=True)
parser.add_argument('--outputdir', help='path to output dir',required=True)
args = parser.parse_args()

chromosome=args.chromosome
haplotype1_read_dictfile=args.haplotype1_read_dict
haplotype2_read_dictfile=args.haplotype2_read_dict

haplotype1_hic_snp_phased_readsfile=args.haplotype1_hic_snp_phased_reads
haplotype2_hic_snp_phased_readsfile=args.haplotype2_hic_snp_phased_reads

haplotype1_contig_lenghtsfile=args.haplotype1_contig_lenghts
haplotype2_contig_lenghtsfile=args.haplotype2_contig_lenghts

outputdir=args.outputdir


def weighted_choice(objects, weights):
    """ returns randomly an element from the sequence of 'objects', 
        the likelihood of the objects is weighted according 
        to the sequence of 'weights', i.e. percentages."""
    weights = np.array(weights, dtype=np.float64)
    sum_of_weights = weights.sum()
    # standardization:
    np.multiply(weights, 1 / sum_of_weights, weights)
    weights = weights.cumsum()
    x = random()
    for i in range(len(weights)):
        if x < weights[i]:
            return objects[i]


#save reads
def save_list(MyList, output_file_name):
        MyFile=open(output_file_name,'w')
        MyList=map(lambda x:x+'\n', MyList)
        MyFile.writelines(MyList)
        MyFile.close()



haplotype2_read_dict = np.load(haplotype2_read_dictfile).item()
haplotype1_read_dict = np.load(haplotype1_read_dictfile).item()

#load in phased reads 
haplotype2_hic_reads = list(pd.read_csv(haplotype2_hic_snp_phased_readsfile, sep=' ', header=None)[0])
haplotype1_hic_reads = list(pd.read_csv(haplotype1_hic_snp_phased_readsfile, sep=' ', header=None)[0])
#turn phased reads into a dictionaty - its faster
haplotype2_hic_reads_dict={}
for read in haplotype2_hic_reads:
        haplotype2_hic_reads_dict[read]=1

#turn phased reads into a dictionaty - its faster
haplotype1_hic_reads_dict={}
for read in haplotype1_hic_reads:
        haplotype1_hic_reads_dict[read]=1

haplotype2_input  = list(haplotype2_read_dict.keys())
haplotype2_size = int(round(len(haplotype2_input)/19))
haplotype2_chunks = [haplotype2_input[i:i+int(haplotype2_size)] for i in range(0, len(haplotype2_input), int(haplotype2_size))]  

average_insert_size = 967.6

#make a list of contigs
haplotype2_contig_list = list(pd.read_csv(haplotype2_contig_lenghtsfile, sep='\t', header=None)[0])
haplotype1_contig_list = list(pd.read_csv(haplotype1_contig_lenghtsfile, sep='\t', header=None)[0])

#make a list of contig lengths 
haplotype1_contig_len = pd.read_csv(haplotype1_contig_lenghtsfile, sep='\t', header=None)
haplotype1_contig_len_dict={}
for item in range(len(list(haplotype1_contig_len[0]))):
	haplotype1_contig_len_dict[haplotype1_contig_len[0][item]]=int(haplotype1_contig_len[1][item])


haplotype2_contig_len = pd.read_csv(haplotype2_contig_lenghtsfile, sep='\t', header=None)
haplotype2_contig_len_dict={}
for item in range(len(list(haplotype2_contig_len[0]))):
	haplotype2_contig_len_dict[haplotype2_contig_len[0][item]]=int(haplotype2_contig_len[1][item])

#make the probabilities of each catagory 
distance_groups = [1000, 1000000, 3162278, 5623413,1000000000]
weights = [0.54, 0.34, 0.044, 0.014, 0.064]
distance_groups_dict = {}
for item in range(len(distance_groups)):
    distance_groups_dict[distance_groups[item]]=weights[item]

poor_quality_reads = []
equal_mapping_reads=[]
equal_mapping_reads_pos=[]
haplotype1_reads=[]
haplotype2_reads=[]
def split_reads_parallel(haplotype2_chunks):
        chunk = haplotype2_chunks
        print("sent to core")
        poor_quality_reads = []
        equal_mapping_reads_pos=[]
        equal_mapping_reads=[]
        haplotype1_reads=[]
        haplotype2_reads=[]
        for haplotype2_read in chunk:
                if haplotype2_read in haplotype1_read_dict:
                #if the hetreozygous snps indicate that the reads should be in a specific group:
                        if haplotype2_read in haplotype2_hic_reads_dict and haplotype2_read not in haplotype1_hic_reads_dict:
                                haplotype2_reads.append(haplotype2_read)
                        elif haplotype2_read in haplotype1_hic_reads_dict and haplotype2_read not in haplotype2_hic_reads_dict:
                                haplotype1_reads.append(haplotype2_read)
                        #if all reads are mapped
                        elif haplotype2_read_dict[haplotype2_read]['mapping_read']=='mapped' and haplotype2_read_dict[haplotype2_read]['mapping_mate']=='mapped' and haplotype1_read_dict[haplotype2_read]['mapping_read']=='mapped' and haplotype1_read_dict[haplotype2_read]['mapping_mate']=='mapped':
                                #if all the mapping is bad 
                                if int(haplotype2_read_dict[haplotype2_read]['read_mapq'])<10 and int(haplotype2_read_dict[haplotype2_read]['mate_mapq'])<10 and int(haplotype1_read_dict[haplotype2_read]['read_mapq'])<10 and int(haplotype1_read_dict[haplotype2_read]['mate_mapq'])<10:
                                        poor_quality_reads.append(haplotype2_read)
                                #if the haplotype1 mapping is bad but the haplotype2 is not 
                                elif int(haplotype2_read_dict[haplotype2_read]['read_mapq'])>10 and int(haplotype2_read_dict[haplotype2_read]['mate_mapq'])>10 and int(haplotype1_read_dict[haplotype2_read]['read_mapq'])<10 and int(haplotype1_read_dict[haplotype2_read]['mate_mapq'])<10:
                                        haplotype2_reads.append(haplotype2_read)
                                #if the haplotype2 mapping is bad but the wild type is not 
                                elif int(haplotype2_read_dict[haplotype2_read]['read_mapq'])<10 and int(haplotype2_read_dict[haplotype2_read]['mate_mapq'])<10 and int(haplotype1_read_dict[haplotype2_read]['read_mapq'])>10 and int(haplotype1_read_dict[haplotype2_read]['mate_mapq'])>10:
                                        haplotype1_reads.append(haplotype2_read)
                                #if the haplotype2 read mapping is better than the haplotype1 read mapping and none of the haplotype2 mapping is 0
                                elif int(haplotype2_read_dict[haplotype2_read]['read_mapq']) + int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) > int(haplotype1_read_dict[haplotype2_read]['read_mapq']) + int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):
                                        haplotype2_reads.append(haplotype2_read)
                                #if the haplotype1 read mapping is better than the haplotype2 read mapping and none of the haplotype2 mapping is 0
                                elif int(haplotype2_read_dict[haplotype2_read]['read_mapq']) + int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) < int(haplotype1_read_dict[haplotype2_read]['read_mapq']) + int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):
                                        haplotype1_reads.append(haplotype2_read)
                                #if the mapping is the same 
                                elif int(haplotype2_read_dict[haplotype2_read]['read_mapq']) + int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) == int(haplotype1_read_dict[haplotype2_read]['read_mapq']) + int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):
                                #if reads map to the same contig
                                        if haplotype2_read_dict[haplotype2_read]['read_contig']==haplotype2_read_dict[haplotype2_read]['mate_contig'] and haplotype1_read_dict[haplotype2_read]['read_contig']==haplotype1_read_dict[haplotype2_read]['mate_contig']:
                                                haplotype2_insert_size = abs(int(haplotype2_read_dict[haplotype2_read]['read_pos'])-int(haplotype2_read_dict[haplotype2_read]['mate_pos']))
                                                haplotype1_insert_size = abs(int(haplotype1_read_dict[haplotype2_read]['read_pos'])-int(haplotype1_read_dict[haplotype2_read]['mate_pos']))
                                                if haplotype2_insert_size > 7500000 and haplotype1_insert_size < 7500000:
                                                        haplotype1_reads.append(haplotype2_read)
                                                elif haplotype2_insert_size < 7500000 and haplotype1_insert_size > 7500000:
                                                        haplotype2_reads.append(haplotype2_read)
                                                elif haplotype2_insert_size > 7500000 and haplotype1_insert_size > 7500000:
                                                        equal_mapping_reads.append(haplotype2_read)
                                                        equal_mapping_reads_pos.append(haplotype2_read_dict[haplotype2_read]['read_pos'])
                                                else:
                                                        haplotype2_insert_size_rounded = min([ i for i in distance_groups if i >= haplotype2_insert_size], key=lambda x:abs(x-haplotype2_insert_size))
                                                        haplotype1_insert_size_rounded = min([ i for i in distance_groups if i >= haplotype1_insert_size], key=lambda x:abs(x-haplotype1_insert_size))
                                                        if haplotype2_insert_size_rounded==haplotype1_insert_size_rounded:
                                                                equal_mapping_reads.append(haplotype2_read)
                                                                equal_mapping_reads_pos.append(haplotype2_read_dict[haplotype2_read]['read_pos'])
                                                        else:
                                                                assignment = weighted_choice([haplotype2_insert_size_rounded,haplotype1_insert_size_rounded], [distance_groups_dict[haplotype2_insert_size_rounded],distance_groups_dict[haplotype1_insert_size_rounded]])
                                                                if assignment==haplotype2_insert_size_rounded:
                                                                        haplotype2_reads.append(haplotype2_read)
                                                                else:
                                                                        haplotype1_reads.append(haplotype2_read)
                                        #if the chrothrip reads map to the same contig and the wild type do not but the haplotype1 reads are close to ends of contigs
                                        elif haplotype2_read_dict[haplotype2_read]['read_contig']==haplotype2_read_dict[haplotype2_read]['mate_contig'] and haplotype1_read_dict[haplotype2_read]['read_contig']!=haplotype1_read_dict[haplotype2_read]['mate_contig'] and \
                                        (haplotype1_contig_len_dict[haplotype1_read_dict[haplotype2_read]['read_contig']]-int(haplotype1_read_dict[haplotype2_read]['read_pos'])<average_insert_size or int(haplotype1_read_dict[haplotype2_read]['read_pos'])<average_insert_size) and (haplotype1_contig_len_dict[haplotype1_read_dict[haplotype2_read]['mate_contig']]-int(haplotype1_read_dict[haplotype2_read]['mate_pos'])<average_insert_size or int(haplotype1_read_dict[haplotype2_read]['mate_pos'])<average_insert_size):
                                                haplotype1_reads.append(haplotype2_read)
                                        #if the haplotype1 reads map to the same contig and the haplotype2 do not but the haplotype2 reads are close to ends of contigs
                                        elif haplotype2_read_dict[haplotype2_read]['read_contig']!=haplotype2_read_dict[haplotype2_read]['mate_contig'] and haplotype1_read_dict[haplotype2_read]['read_contig']==haplotype1_read_dict[haplotype2_read]['mate_contig'] and \
                                        (haplotype2_contig_len_dict[haplotype2_read_dict[haplotype2_read]['read_contig']]-int(haplotype2_read_dict[haplotype2_read]['read_pos'])<average_insert_size or int(haplotype2_read_dict[haplotype2_read]['read_pos'])<average_insert_size) and (haplotype2_contig_len_dict[haplotype2_read_dict[haplotype2_read]['mate_contig']]-int(haplotype2_read_dict[haplotype2_read]['mate_pos'])<average_insert_size or int(haplotype2_read_dict[haplotype2_read]['mate_pos'])<average_insert_size):
                                              haplotype2_reads.append(haplotype2_read)
                                        #if all reads map to different contigs
                                        else:
                                                if (haplotype2_contig_len_dict[haplotype2_read_dict[haplotype2_read]['read_contig']]-int(haplotype2_read_dict[haplotype2_read]['read_pos'])<average_insert_size or int(haplotype2_read_dict[haplotype2_read]['read_pos'])<average_insert_size) and (haplotype2_contig_len_dict[haplotype2_read_dict[haplotype2_read]['mate_contig']]-int(haplotype2_read_dict[haplotype2_read]['mate_pos'])<average_insert_size or int(haplotype2_read_dict[haplotype2_read]['mate_pos'])<average_insert_size):
                                                        haplotype2_reads.append(haplotype2_read)
                                                elif (haplotype1_contig_len_dict[haplotype1_read_dict[haplotype2_read]['read_contig']]-int(haplotype1_read_dict[haplotype2_read]['read_pos'])<average_insert_size or int(haplotype1_read_dict[haplotype2_read]['read_pos'])<average_insert_size) and (haplotype1_contig_len_dict[haplotype1_read_dict[haplotype2_read]['mate_contig']]-int(haplotype1_read_dict[haplotype2_read]['mate_pos'])<average_insert_size or int(haplotype1_read_dict[haplotype2_read]['mate_pos'])<average_insert_size):
                                                        haplotype1_reads.append(haplotype2_read)
                                                else:
                                                        equal_mapping_reads.append(haplotype2_read)
                                                        equal_mapping_reads_pos.append(haplotype2_read_dict[haplotype2_read]['read_pos'])
                                else:
                                        print('you missed something in mapq', haplotype2_read_dict[haplotype2_read]['read_mapq'], haplotype2_read_dict[haplotype2_read]['mate_mapq'], haplotype1_read_dict[haplotype2_read]['read_mapq'], haplotype1_read_dict[haplotype2_read]['mate_mapq'])
                        #if only the haplotype2 reads map - dont worry about mapq
                        elif (haplotype2_read_dict[haplotype2_read]['mapping_read']=='mapped' or haplotype2_read_dict[haplotype2_read]['mapping_mate']=='mapped') and (haplotype1_read_dict[haplotype2_read]['mapping_read']=='unmapped') and (haplotype1_read_dict[haplotype2_read]['mapping_mate']=='unmapped'):
                                haplotype2_reads.append(haplotype2_read)
                        #if only the haplotype1_reads map - dont worry about mapq
                        elif (haplotype1_read_dict[haplotype2_read]['mapping_read']=='mapped' or haplotype1_read_dict[haplotype2_read]['mapping_mate']=='mapped') and (haplotype2_read_dict[haplotype2_read]['mapping_read']=='unmapped') and (haplotype2_read_dict[haplotype2_read]['mapping_mate']=='unmapped'):
                                haplotype1_reads.append(haplotype2_read)
                        #if only one of the read pair maps in both haplotype2 and haplotype1
                        elif haplotype2_read_dict[haplotype2_read]['mapping_read']=='mapped' and haplotype1_read_dict[haplotype2_read]['mapping_read']=='mapped' and  haplotype2_read_dict[haplotype2_read]['mapping_mate']=='unmapped' and haplotype1_read_dict[haplotype2_read]['mapping_mate']=='unmapped':
                                if int(haplotype2_read_dict[haplotype2_read]['read_mapq']) > int(haplotype1_read_dict[haplotype2_read]['read_mapq']):
                                        haplotype2_reads.append(haplotype2_read)             
                                if int(haplotype2_read_dict[haplotype2_read]['read_mapq']) < int(haplotype1_read_dict[haplotype2_read]['read_mapq']):
                                        haplotype1_reads.append(haplotype2_read)
                                if int(haplotype2_read_dict[haplotype2_read]['read_mapq']) == int(haplotype1_read_dict[haplotype2_read]['read_mapq']):
                                        equal_mapping_reads.append(haplotype2_read)
                                        equal_mapping_reads_pos.append(haplotype2_read_dict[haplotype2_read]['read_pos'])
                        #if only one of the read pair maps in both haplotype2 and haplotype1
                        elif (haplotype2_read_dict[haplotype2_read]['mapping_mate']=='mapped' or haplotype1_read_dict[haplotype2_read]['mapping_mate']=='mapped') and haplotype2_read_dict[haplotype2_read]['mapping_read']=='unmapped' and haplotype1_read_dict[haplotype2_read]['mapping_read']=='unmapped':
                                if int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) > int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):
                                        haplotype2_reads.append(haplotype2_read)             
                                if int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) < int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):
                                        haplotype1_reads.append(haplotype2_read)
                                if int(haplotype2_read_dict[haplotype2_read]['mate_mapq']) == int(haplotype1_read_dict[haplotype2_read]['mate_mapq']):     
                                        equal_mapping_reads.append(haplotype2_read)
                                        equal_mapping_reads_pos.append(haplotype2_read_dict[haplotype2_read]['read_pos'])
                        #if all reads but one map 
                        elif (haplotype2_read_dict[haplotype2_read]['mapping_read']=='mapped' and haplotype2_read_dict[haplotype2_read]['mapping_mate']=='mapped') and (haplotype1_read_dict[haplotype2_read]['mapping_read']=='unmapped' or haplotype1_read_dict[haplotype2_read]['mapping_mate']=='unmapped'):
                                haplotype2_reads.append(haplotype2_read)
                        #if only the haplotype1_reads map - dont worry about mapq
                        elif (haplotype1_read_dict[haplotype2_read]['mapping_read']=='mapped' and haplotype1_read_dict[haplotype2_read]['mapping_mate']=='mapped') and (haplotype2_read_dict[haplotype2_read]['mapping_read']=='unmapped') or (haplotype2_read_dict[haplotype2_read]['mapping_mate']=='unmapped'):
                                haplotype1_reads.append(haplotype2_read)
                        #no mapping for any of the read
                        elif (haplotype2_read_dict[haplotype2_read]['mapping_read']=='unmapped' or haplotype2_read_dict[haplotype2_read]['mapping_mate']=='unmapped') and haplotype1_read_dict[haplotype2_read]['mapping_read']=='unmapped' and haplotype1_read_dict[haplotype2_read]['mapping_mate']=='unmapped':
                                poor_quality_reads.append(haplotype2_read)
                        else:
                                poor_quality_reads.append(haplotype2_read)
                else:
                        #read missing in haplotype1 dict
                        haplotype2_reads.append(haplotype2_read)
        
        return poor_quality_reads, haplotype2_reads, haplotype1_reads, equal_mapping_reads, equal_mapping_reads_pos



pool = mp.Pool(9) #change if you change the number of cores the job is submitted to
output = pool.map(split_reads_parallel, haplotype2_chunks)
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


for haplotype1_read in haplotype1_read_dict:
        if haplotype1_read not in haplotype2_read_dict:
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
save_list(haplotype1_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype1_reads_all_',chromosome,'_',sample,'.txt'))
save_list(haplotype2_reads, "{}{}{}{}{}{}".format(outputdir,'haplotype2_reads_all_',chromosome,'_',sample,'.txt'))

