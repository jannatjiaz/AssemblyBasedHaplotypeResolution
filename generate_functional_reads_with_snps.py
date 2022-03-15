import pandas as pd 
import numpy as np
import multiprocessing as mp
import os 
import sys
from os import system
import re
import csv 
import argparse 

parser = argparse.ArgumentParser(description='usage:  generate_functional_reads_with_snps.py --chromosome chromosome --run run --haplotype1_snps haplotype1_snps --haplotype2_snps haplotype2_snps --mpileup_reads mpileup_reads --outputdir path')
parser.add_argument('--chromosome', help='chromosome',required=True)
parser.add_argument('--run', help='hic run number',required=True)
parser.add_argument('--haplotype1_snps', help='snps from haplotype 1 generated from e.g. ccs',required=True)
parser.add_argument('--haplotype2_snps', help='snps from haplotype 2 generated from e.g. ccs',required=True)
parser.add_argument('--mpileup_reads', help='path to mpileup_with_reads_hic.py output',required=True)
parser.add_argument('--outputdir', help='path to the output directory',required=True)
args = parser.parse_args()


chromosome=args.chromosome
run=args.run
haplotype1_snps_file=args.haplotype1_snps
haplotype2_snps_file=args.haplotype2_snps
mpileup_reads_file=args.mpileup_reads
outputdir=args.outputdir

#save reads
def save_list(MyList, output_file_name):
        MyFile=open(output_file_name,'w')
        MyList=map(lambda x:x+'\n', MyList)
        MyFile.writelines(MyList)
        MyFile.close()

#load snp list
haplotype2_snps=list(pd.read_csv(haplotype2_snps_file,sep="\t", header=None)[0])
haplotype2_snps_dict={x:1 for x in haplotype2_snps} 

haplotype1_snps=list(pd.read_csv(haplotype1_snps_file,sep="\t", header=None)[0])
haplotype1_snps_dict={x:1 for x in haplotype1_snps} 

#load mpile up with reads file
mpile_up_input_hic=pd.read_csv(mpileup_reads_file,sep="\t",header=None)
mpile_up_input_hic.columns=['chr','pos','ref_snp','read','read_query_pos','actual_snp']

snp_id = list(mpile_up_input_hic['chr'].astype(str)+'_'+mpile_up_input_hic['pos'].astype(str)+'_'+mpile_up_input_hic['actual_snp'].astype(str)+"_"+mpile_up_input_hic['read'].astype(str))
#make the mpileup dictionary where each key is an SV id and each key contains, chr, pos, actual_snp, read
mpile_up_input_hic_dict={x:{} for x in snp_id} #make an empty dictionary for each of the SVs
for item in range(len(snp_id)):
        mpile_up_input_hic_dict[snp_id[item]]= {"chr":mpile_up_input_hic['chr'][item], "pos":mpile_up_input_hic["pos"][item], "ref_snp":mpile_up_input_hic["ref_snp"][item],"read":mpile_up_input_hic["read"][item],"read_query_pos":mpile_up_input_hic["read_query_pos"][item], "actual_snp":mpile_up_input_hic["actual_snp"][item]}

print(len(mpile_up_input_hic_dict))


haplotype2_reads={}
haplotype1_reads={}
#iterate over the snp_ids
for snp in mpile_up_input_hic_dict:
    #if there is a base at that posision, ie its not a deltion
    if mpile_up_input_hic_dict[snp]['read_query_pos'] != 'None':
        #generate the speicific snp id
        snp_identity=str(mpile_up_input_hic_dict[snp]['chr'])+'_'+str(mpile_up_input_hic_dict[snp]['pos'])+'_'+str(mpile_up_input_hic_dict[snp]['actual_snp'])
        #if snp is in haplotype1itpic snps:
        if snp_identity in haplotype1_snps_dict:
            #add the snp id to the haplotype1 snps list if its not already there 
            if snp_identity not in haplotype1_reads:
                haplotype1_reads[mpile_up_input_hic_dict[snp]['read']]=1
       #if snp is in haplotype2 snps:
        if snp_identity in haplotype2_snps_dict:
            #add the snp id to the haplotype2 snps list if its not already there 
            if snp_identity not in haplotype2_reads:
                haplotype2_reads[mpile_up_input_hic_dict[snp]['read']]=1

print(len(haplotype2_reads))
print(len(haplotype1_reads))

##make sure all snps are unambigiously phased
haplotype2_reads_final=[]
for key in haplotype2_reads:
    if key not in haplotype1_reads:
        haplotype2_reads_final.append(key)

haplotype1_reads_final=[]
for key in haplotype1_reads:
    if key not in haplotype2_reads:
        haplotype1_reads_final.append(key)



print(len(haplotype1_reads_final))
print(len(haplotype2_reads_final))

save_list(haplotype2_reads_final, "{}{}{}{}{}{}".format(outputdir,'/haplotype2_reads_',run,'_',chromosome,'.txt'))
save_list(haplotype1_reads_final, "{}{}{}{}{}{}".format(outputdir,'/haplotype1_reads_',run,'_',chromosome,'.txt'))
