#code based on: http://pysam.readthedocs.org/en/latest/
#argparse info: http://www.cyberciti.biz/faq/python-command-line-arguments-argv-example/
import pysam
import argparse
parser = argparse.ArgumentParser(description='usage:  samtools_view.py  --bam reads.bam --fasta reference_genomic.fasta  --chr ctg2_right --start 100140')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--fasta', help='Input reference fasta file name',required=True)
parser.add_argument('--chr',help='target contig or chromosome', required=True)
parser.add_argument('--start',help='target region start (bp)', required=True)
args = parser.parse_args()
args.start = int(args.start) -1  #pysam is 0-based, but human minds are not
samfile = pysam.AlignmentFile(args.bam, "rb")
fastafile = pysam.FastaFile( args.fasta )
for pileupcolumn in samfile.pileup(args.chr, args.start , args.start+1, min_base_quality=0):  #pysam is 0-based, but human minds are not
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del and pileupread.indel==0 and (pileupcolumn.pos >= args.start ) and (pileupcolumn.pos <= args.start):
            print ('%s\t%s\t%s\t%s\t%s\t%s' %
                 (args.chr,
                  pileupcolumn.pos +1,                                                     # +1 needed because pysam is 0-based and most programs (and human heads) are 1-based
                  fastafile.fetch(args.chr,pileupcolumn.reference_pos ,pileupcolumn.reference_pos +1), 
                  pileupread.alignment.query_name,
                  pileupread.query_position,  
                  "*"))     
        if not pileupread.is_del and pileupread.indel==0 and (pileupcolumn.pos >= args.start ) and (pileupcolumn.pos <= args.start  ):
           print ('%s\t%s\t%s\t%s\t%s\t%s' %
                 (args.chr,
                  pileupcolumn.pos +1,
                  fastafile.fetch(args.chr,pileupcolumn.reference_pos ,pileupcolumn.reference_pos +1),  
                  pileupread.alignment.query_name,
                  pileupread.query_position,  
                  pileupread.alignment.query_sequence[pileupread.query_position]))
        if not pileupread.is_del and pileupread.indel>0 and (pileupcolumn.pos >= args.start ) and (pileupcolumn.pos <=args.start ):
             print ('%s\t%s\t%s\t%s\t%s\t%s' %
                 (args.chr,
                  pileupcolumn.pos +1,
                  fastafile.fetch(args.chr,pileupcolumn.reference_pos ,pileupcolumn.reference_pos +1),  
                  pileupread.alignment.query_name,
                  pileupread.query_position,  
                  pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1] )) #+1 because python goes up to the last number in a list
        if not pileupread.is_del and pileupread.indel<0 and (pileupcolumn.pos >= args.start ) and (pileupcolumn.pos <=args.start ):
            print ('%s\t%s\t%s\t%s\t%s\t%s' %
                 (args.chr,
                  pileupcolumn.pos +1,
                  fastafile.fetch(args.chr,pileupcolumn.reference_pos ,pileupcolumn.reference_pos +1),  
                  pileupread.alignment.query_name,
                  pileupread.query_position,  
                  fastafile.fetch(args.chr,pileupcolumn.reference_pos ,pileupcolumn.reference_pos-pileupread.indel +1))) #+1 because python goes up to the last number in a list
samfile.close()
fastafile.close()
