# HapFunc
Haplotype resolve epigenetic and transcriptomic data types

HapFunc is used to haplotpye resolve epigentic and transcriptomic datasets based on haplotpye specific assemblies. Each code is tailored to features of the underlying dataset. Resolution can be done on paired-end illumina sequencing (e.g. ChIP-seq or ATAC-seq), paired-end illumina sequencing from Hi-C or long reads (e.g. iso-seq).

## Generate SNPs in CCS reads

In order to run this code you need a SNP file for haplotype A and haplotype B derived from CCS reads. You need a list of the position of variants that were used for whatshap phasing (filter to remove simple repeats regions). You can split this list into multiple files to parallelise. Then run mpileup_CCS.py:

> while read VARIANT; do \
>       python mpileup_CCS.py --bam CCS.sorted.bam --fasta genome.fa  --chr chromosome --start $VARIANT  >> mpileup_CCS.txt \
> done < positions.txt 

The output of this will be: chromosome, position, ref_base, read_name, postion_in_read, base

Once you have the CCS list you can determine which SNPs belong to each haplotype. You will need a list of all the CCS reads that are from haplotype 1 and all the reads from haplotype 2 (haplotype1_readnames.txt and haplotype2_readnames.txt)

> python generate_snps_from_ccs.py  \
> --chromosome chromosome
> --haplotype1_reads haplotype1_readnames.txt \
> --haplotype2_reads haplotype2_readnames.txt \
> --mpileup_reads mpileup_CCS.txt \
> --outputdir output_directory

This will produce two files e.g. haploype1_snps_chr1.txt and haploype2_snps_chr1.txt. These will contain the snps at each position that has been inputted for the chromothrpitic and wild-type alleles.

## Haplotype resolve epigenetic and transcriptomic reads

Now the SNPs for each haplotpye have been determined, you will be able to phase the epigenetic and transcriptomic datasets. The same steps are needed regardless of data type but the underlying rules change so each script has been tailored for the specific datasets. Here instructions are written for haplotpye resolving Hi-C reads however the scripts for haplotype resolving other functional dataset can be replaces and the parameters stay the same.  

### Generate SNPs in epigenetic and transcriptomic reads

The first step is to determine which SNPs are in which reads, similar to the process used for CCS reads. Use the same position file as was used for determining SNPs in CCS reads. Run mpileup_hic.py

> while read VARIANT; do 
>   python mpileup_hic.py  --bam hic.sorted.bam --fasta genome.fa --chr chromosome --start $VARIANT  >> mpileup_hic.txt
> done <  positions.txt

Now you have the which SNPs are present in which reads, you can determine which reads are from which haplotpye based on SNPs. Run generate_hic_reads_with_snps.py. For example for chr1:

> python generate_hic_reads_with_snps.py \
> --chromosome chr1 --hic_run run1 \
> --haplotype1_snps haplotype1_snps_chr1.txt \
> --haplotype2_snps haplotype2_snps_chr2.txt \
> --mpileup_reads mpileup_hic.txt \
> --outputdir output_directory

### Haplotype resolution 

The majority of haplotype resolution is from SNP presence but some resolution can be gained from differences in mapping. For Hi-C reads, further information is gained from a weighted probability that two reads are interacting at specific distances. To take these factors into account, you will need a list of the readnames for all the 



