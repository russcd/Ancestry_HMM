##This is a simple updated version of our vcf2ahmm script to take an input vcf file and output a file compatible ancestry_hmm. 

The script takes the following required options.

python3 vcf2ahmm.py -v [vcf_file] -s [sample2population mappings] > [ahmm_input_file] 

Where the sample2population mapping file is a text-based table with two columns. 

1. sample id exactly as it appears in the vcf file
2. the population to which that individual belongs. This can be either one of the ancestral populations indicated with an integer (0,1,2..k), or the sample is admixed in which case the population must read "admixed"

e.g. 

sample1	0
sample2	0
sample3 1
sample4 1
sample5 admixed

This file has two individuals from ancestral population 0, two individuals from ancestral population 1, and a single admixed sample.

##Assumptions:

1. Sample are either diploid or haploid. The allele depth field (AD), immediately follows the genotype (GT) in the vcf format. This is typical of most VCF files. 
2. A single uniform per basepair, per generation recombination rate should be used. Default: 1e-8. 

##Optional arguments:

1. "-g 1" indicates that admixed sample genotypes (rather than allele counts) should be used
2. "-r [float]" sets the 
3. "-m [int]" minimum distance in basepairs between successful SNPs to be included in this analysis 
4. "--min_total [int]" minimum number of samples in each ancestral population to consider a site. Default 10. 
5. "--min_diff [float]" minimum allele frequency difference between any pair of ancestral populations to include a site. I.e., this selects AIMs. Default 0.1. 
