

import pandas as pd
import pysam 
import argparse

parser = argparse.ArgumentParser(description='usage:  generate_haplotype_sam.py  \
    --chromosome chromosome \
    --sample_name sample_name\
    --haplotype haplotype\
    --outputdir path\
    --inputdir path\
    --bam bamfile')
parser.add_argument('--chromosome', help='chromosome name',required=True)
parser.add_argument('--haplotype', help='1 or 2 depending on the haplotype',required=True)
parser.add_argument('--sample_name', help='name of sample',required=True)
parser.add_argument('--outputdir', help='path to output dir',required=True)
parser.add_argument('--inputdir', help='path to input reads text file  dir',required=True)
parser.add_argument('--bam', help='path to all aligned bam reads ',required=True)
args = parser.parse_args()

chromosome=args.chromosome
sample_name=args.sample_name
haplotype=args.haplotype
outputdir=args.outputdir
inputdir=args.inputdir
bamfile=args.bam

bam = pysam.AlignmentFile(bamfile, "rb") #load bam file

obam = pysam.AlignmentFile("{}{}{}{}{}{}{}{}".format(outputdir,"/haplotype",haplotype,"_reads_assigned_",chromosome,"_",sample_name,".sam"), "w", template=bam)
haplotype_read_names = pd.read_csv("{}{}{}{}{}{}{}{}".format(inputdir,"/haplotype",haplotype,"_reads_assigned_",chromosome,"_",sample_name,".txt"))
haplotype_read_names.columns =['readnames']
haplotype_read_names= list(haplotype_read_names['readnames'])
haplotype_read_dict={}
for read in haplotype_read_names:
        haplotype_read_dict[read]=1

for b in bam.fetch():
    if b.query_name in haplotype_read_dict:
        obam.write(b)

obam.close()


outbam = pysam.AlignmentFile("{}{}{}{}{}{}{}{}".format(outputdir,"/haplotype",haplotype,"_reads_all_",chromosome,"_",sample_name,".sam"), "w", template=bam)
haplotype_read_names = pd.read_csv("{}{}{}{}{}{}{}{}".format(inputdir,"/haplotype",haplotype,"_reads_",chromosome,"_",sample_name,".txt"))
haplotype_read_names.columns =['readnames']
haplotype_read_names= list(haplotype_read_names['readnames'])
haplotype_read_dict={}
for read in haplotype_read_names:
        haplotype_read_dict[read]=1

for b in bam.fetch():
    if b.query_name in haplotype_read_dict:
        outbam.write(b)

outbam.close()

bam.close()



