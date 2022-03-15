import numpy as np
import argparse 

parser = argparse.ArgumentParser(description='usage:  create_haplotype1_dict.py  \
    --haplotype_reads haplotype_reads \
    --haplotype 1 \
    --outputdir path')
parser.add_argument('--haplotype_reads', help='hap1 or hap2 alignments from previous step',required=True)
parser.add_argument('--haplotype', help='1 or 2 depending on the haplotype',required=True)
parser.add_argument('--outputdir', help='path to output dir',required=True)
args = parser.parse_args()

haplotype_readsfile=args.haplotype_reads
haplotype_number=args.haplotype
outputdir=args.outputdir


keys = ["read","mapping_read","mapping_mate","read_mapq","mate_mapq","read_contig","mate_contig","read_pos","mate_pos"]
#file = open('../final_split/hic1/chr20/hic1_chr20_hap2_alignments.txt', 'r')
file = open(haplotype_readsfile, 'r')
results = []
while True:
    line = file.readline()
    if not line:
        break
    else:
        content = line.rstrip("\n").split(',')
        dict = {}
        index = 0
        for value in content:
            if value != '\n':
                pair = {keys[index]: value}
                dict.update(pair)
                index += 1
        if dict != {}:                          # Prevent empty lines from appending to results
            results.append(dict)

new_dict = {item['read']:item for item in results}

np.save("{}{}".format(outputdir,'/haplotype',haplotype_number,'_read_dict.npy'), new_dict) 