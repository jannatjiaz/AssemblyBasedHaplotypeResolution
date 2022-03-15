import pandas as pd
import numpy as np
import argparse 

parser = argparse.ArgumentParser(description='usage:  generate_ccs_snps.py --chromosome chromosome --haplotype1_reads haplotype1_reads --haplotype2_reads haplotype2_reads --mpileup_reads mpileup_reads --outputdir path')
parser.add_argument('--chromosome', help='chromosome',required=True)
parser.add_argument('--haplotype1_reads', help='readnames haplotype1 ccs',required=True)
parser.add_argument('--haplotype2_reads', help='readnames haplotype2 ccs',required=True)
parser.add_argument('--mpileup_reads', help='path mpileup_with_reads.py output',required=True)
parser.add_argument('--outputdir', help='path to the output directory',required=True)
args = parser.parse_args()


chromosome=args.chromosome
haplotype1_reads_file=args.haplotype1_reads
haplotype2_reads_file=args.haplotype2_reads
mpileup_reads_file=args.mpileup_reads
outputdir=args.outputdir


#save reads
def save_list(MyList, output_file_name):
        MyFile=open(output_file_name,'w')
        MyList=map(lambda x:x+'\n', MyList)
        MyFile.writelines(MyList)
        MyFile.close()

#import the haplotype1 reads and convert it into a list
haplotype1_reads=list(pd.read_csv(haplotype1_reads_file,sep="\t", header=None)[0])
haplotype1_reads_dict={x:1 for x in haplotype1_reads} 
#import the wild type reads and convert it into a list
haplotype2_reads=list(pd.read_csv(haplotype2_reads_file,sep="\t", header=None)[0])
haplotype2_reads_dict={x:1 for x in haplotype2_reads} 

#import mpile_up
mpile_up_input_ccs=pd.read_csv(mpileup_reads_file,sep="\t",header=None)
mpile_up_input_ccs.columns=['chr','pos','ref_snp','read','read_query_pos','actual_snp']

#make the SV dictionary where each key is an SV id and each key contains, SV_id, chr1, pos1, strand1, chr2, pos2, strand2 and reads
snp_id = list(mpile_up_input_ccs['chr'].astype(str)+'_'+mpile_up_input_ccs['pos'].astype(str)+'_'+mpile_up_input_ccs['actual_snp'].astype(str)+"_"+mpile_up_input_ccs['read'].astype(str))
mpile_up_input_ccs_dict={x:{} for x in snp_id} #make an empty dictionary for each of the SVs
for item in range(len(snp_id)):
        mpile_up_input_ccs_dict[snp_id[item]]= {"chr":mpile_up_input_ccs['chr'][item], "pos":mpile_up_input_ccs["pos"][item], "ref_snp":mpile_up_input_ccs["ref_snp"][item],"read":mpile_up_input_ccs["read"][item],"read_query_pos":mpile_up_input_ccs["read_query_pos"][item], "actual_snp":mpile_up_input_ccs["actual_snp"][item]}

print(len(mpile_up_input_ccs_dict))

haplotype2_snps={}
haplotype1_snps={}
for snp in mpile_up_input_ccs_dict:
    #as long as there is a base in that positions, ie its not a deletion
    if mpile_up_input_ccs_dict[snp]['read_query_pos'] != 'None':
        #if the read is in haplotype1 reads list
        if mpile_up_input_ccs_dict[snp]['read'] in haplotype1_reads_dict:
            #generate a snp id
            phase = mpile_up_input_ccs_dict[snp]['chr']+'_'+mpile_up_input_ccs_dict[snp]['pos'].astype(str)+'_'+mpile_up_input_ccs_dict[snp]['actual_snp']
            #add the snp id to the haplotype1 snps list if its not already there 
            if phase not in haplotype1_snps:
                haplotype1_snps[phase]=1
        #if the read in in haplotype2 reads list
        if mpile_up_input_ccs_dict[snp]['read'] in haplotype2_reads_dict:
            #generate a snp id
            phase = mpile_up_input_ccs_dict[snp]['chr']+'_'+mpile_up_input_ccs_dict[snp]['pos'].astype(str)+'_'+mpile_up_input_ccs_dict[snp]['actual_snp']
            #add the snp id to the haplotype2 snps list if its not already there 
            if phase not in haplotype2_snps:
                haplotype2_snps[phase]=1
                #haplotype2_snps.append(phase)

##make sure all snps are unambigiously phased
haplotype2_snps_final=[]
for key in haplotype2_snps:
    if key not in haplotype1_snps:
        haplotype2_snps_final.append(key)

haplotype1_snps_final=[]
for key in haplotype1_snps:
    if key not in haplotype2_snps:
        haplotype1_snps_final.append(key)


#save the lists 
save_list(haplotype2_snps_final, "{}{}{}{}".format(outputdir,'/haplotype2_snps_',chromosome,'.txt'))
save_list(haplotype1_snps_final, "{}{}{}{}".format(outputdir,'/haplotype1_snps_',chromosome,'.txt'))
