# HapFunc
Haplotype resolve epigenetic and transcriptomic data types

HapFunc is used to haplotpye resolve epigentic and transcriptomic datasets based on haplotpye specific assemblies. Each code is tailored to features of the underlying dataset. Resolution can be done on paired-end illumina sequencing (e.g. ChIP-seq or ATAC-seq), paired-end illumina sequencing from Hi-C or long reads (e.g. iso-seq).

## Generate SNPs in CCS reads

In order to run this code you need a SNP file for haplotype A and haplotype B derived from CCS reads. You need a list of the position of variants that were used for whatshap phasing (filter to remove simple repeats regions). You can split this list into multiple files to parallelise. Then run mpileup_CCS.py:

> while read VARIANT; do \
>       python mpileup_CCS.py --bam CCS.sorted.bam --fasta genome.fa  --chr chromosome --start $VARIANT  >> mpileup_CCS.txt \
> done < positions.txt 

The output of this will be: chromosome, position, ref_base, read_name, postion_in_read, base

Once you have the this output, you can determine which SNPs belong to each haplotype. You will need a list of all the CCS reads that are from haplotype 1 and all the reads from haplotype 2 (haplotype1_readnames.txt and haplotype2_readnames.txt)

> python generate_snps_from_ccs.py  \
> --chromosome chromosome
> --haplotype1_reads haplotype1_readnames.txt \
> --haplotype2_reads haplotype2_readnames.txt \
> --mpileup_reads mpileup_CCS.txt \
> --outputdir output_directory

This will produce two files e.g. haploype1_snps_chr1.txt and haploype2_snps_chr1.txt. These will contain the snps at each position that has been inputted for haplotype 1 and haplotype 2.

## Haplotype resolve epigenetic and transcriptomic reads

Now the SNPs for each haplotpye have been determined, you will be able to phase the epigenetic and transcriptomic datasets. The same steps are needed regardless of data type but depending on the data type, the script has been tailored for the unique properties. Here instructions are written for haplotpye resolving Hi-C reads however the scripts for haplotype resolving other functional dataset can be replaced and the input parameters are the same.  

### Generate SNPs in epigenetic and transcriptomic reads

The first step is to determine which SNPs are in which reads, similar to the process used for CCS reads. Use the same position file as was used for determining SNPs in CCS reads. Run mpileup_hic.py

> while read VARIANT; do 
>   python mpileup_hic.py  --bam hic.sorted.bam --fasta genome.fa --chr chromosome --start $VARIANT  >> mpileup_hic.txt
> done <  positions.txt

Please note use the mpileup script that matches the data type. For example for ChIP and ATAC use mpileup_pairedend.py or for long read iso-seq use mpileup_longread.py. The same parameters are used for all scripts.


Now you have the which SNPs are present in which reads, you can determine which reads are from which haplotpye based on SNPs. Run generate_hic_reads_with_snps.py. This can be used on all data types. For example for chr1:

> python generate_functional_reads_with_snps.py \
> --chromosome chr1 --hic_run run1 \
> --haplotype1_snps haplotype1_snps_chr1.txt \
> --haplotype2_snps haplotype2_snps_chr2.txt \
> --mpileup_reads mpileup_hic.txt \
> --outputdir output_directory

This will give you two list of reads: one of which is reads which cover SNPs from haplotype 1 and the other that covers haplotype 2.

### Haplotype resolution 

The majority of haplotype resolution is from SNP presence but some resolution can be gained from differences in mapping. For Hi-C reads, further information is gained from a weighted probability that two reads are interacting at specific distances. You can also take these factors into account in phasing. 

Processing bam files in python is computationally slow and intensitive therefore isolating necessary mapping proporties using bash and samtools speeds the process up.

Fist use samtools (https://github.com/samtools) in order to get the position of mapped reads. You will need to map all the Hi-C reads to the haplotype 1 assembly and all the reads to the haplotype 2 assembly:

> samtools view -F 256 -f65  allreads_hap1.sorted.bam | cut -f1,3,4,5  > R1_hap1.txt \
> samtools view -F 256 -f129  allreads_hap1.sorted.bam | cut -f1,3,4,5 > R2_hap1.txt \
> samtools view -F 256 -f65  allreads_hap2.sorted.bam | cut -f1,3,4,5  > R1_hap2.txt \
> samtools view -F 256 -f129  allreads_hap2.sorted.bam | cut -f1,3,4,5 > R2_hap2.txt 

Now these mapping perameters can be combined using:

> awk 'NR==FNR {a[$1]=$0;next}{print $0, ($1 in a ? a[$1]:"NA")}' R2_hap1.txt R1_hap1.txt | awk '{if ($4!=0 && $8!=0) print $1,"mapped","mapped",$4,$8,$2,$6,$3,$7; else if ($4!=0 && $8==0) print $1,"mapped","unmapped",$4,$8,$2,$6,$3,$7; else if ($4==0 && $8!=0) print $1,"unmapped","mapped",$4,$8,$2,$6,$3,$7; }' | tr " " "," > hap1_alignments.txt

> awk 'NR==FNR {a[$1]=$0;next}{print $0, ($1 in a ? a[$1]:"NA")}' R2_hap2.txt R1_hap2.txt | awk '{if ($4!=0 && $8!=0) print $1,"mapped","mapped",$4,$8,$2,$6,$3,$7; else if ($4!=0 && $8==0) print $1,"mapped","unmapped",$4,$8,$2,$6,$3,$7; else if ($4==0 && $8!=0) print $1,"unmapped","mapped",$4,$8,$2,$6,$3,$7; }' | tr " " "," > hap2_alignments.txt

If using long-reads, for example iso-seq, it is not appropriate to use R1 and R2 so instead use:

> samtools view -F 256 allreads_hap1.sorted.bam | cut -f1,3,4,5 | awk '{print $1 "," "mapped" "," $2 "," $3 "," $4 }' > hap1_alignments.txt \
> samtools view -F 256 allreads_hap2.sorted.bam | cut -f1,3,4,5  | awk '{print $1 "," "mapped" "," $2 "," $3 "," $4 }'> hap2_alignments.txt 

This information is then converted into a python dictionary:







