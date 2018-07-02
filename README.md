# What is Introme?

Introme identifies high quality, rare, intronic, damaging genetic variants resulting from whole genome sequencing and following an expected mendelian inheritance pattern.

It takes a VCF file containing the millions of variants found in a typical family trio and applies a number of filters until a putatively damaging list of intronic variants is left. Further prediction and prioritisation is achieved using database annotations and calculation of potential to create cryptic splice sites.

# Developers

Introme was conceived and initially specced as an idea by Dr. Mark Cowley and Dr. Velimir Gayevskiy. Sarah Beecroft carried out initial implementation. Dr. Velimir Gayevskiy reimplemented and added features.

# Requirements

## You need the following software installed

Introme is developed on macOS and runs with the following dependencies which can be mostly found using the Homebrew package manager for macOS. It should work just fine on Linux with the same dependencies installed as it is a Shell script.

* bedtools
* tabix
* vcfanno
* bcftools
* bgzip
* tabix
* MaxEntScan (supplied with Introme from http://genes.mit.edu/burgelab/maxent/download/)

## The VCF file supplied should be

* Created using GATK HaplotypeCaller (other variant callers may work but are currently untested)
* Decomposed and normalised (i.e. no multi-allelic variants) with vt
* Annotated with VEP (optional but very useful for knowing what genes the variants are in)

## Annotations

* CADD v1.3 VCF created using the instructions at: https://github.com/brentp/vcfanno/blob/master/docs/examples/cadd.md
* gnomad.genomes.sites.merged.AF\_AC_AN\_only.vcf.gz
* MGRB variant allele frequencies (these are not public so remove this annotation manually)
* Branchpointer (Machine-learning annotation of human splicing branchpoints) from https://osf.io/hrqvq/
* GENCODE coordinates to gene names
* dbscSNV
* SPIDEX

# Using Introme

## Options

* -a -- one or more affected sample names (specify multiple times if more than one), must have the same sample name as in the VCF file
* -b -- path to search space BED file
* -i -- inheritance pattern, can be one of "denovo" or "autrec"
* -p -- output prefix
* -v -- path to VCF file
* -r -- path to reference genome fasta file (must be indexed with samtools)

## Example

./run\_introme.sh -b subsetting/UCSC\_introns.bed.gz -v input/Fam1_jointcall.hc.vqsr.decomposed.normalised.vep.vcf.gz -p Fam1 -i denovo -a A001C -r annotations/hs37d5.fasta-index/genome.fa

# MaxEntScan consequence logic

MaxEntScan produces a numeric score for a given string of 9 bases for 5' splice sites and 23 bases for 3' splice sites. This score denotes the affinity of the bases for acting as a splice site. However, in the context of determining whether a single SNV in an intron is creating a new splice site, it is not obvious whether the splice site created is a 5' or 3' splice site and also where within the splice site the variant falls.

Introme attempts to find the maximal splice site creation potential of the variant. It does this by first supposing that the variant is the first base of a new splice site and calculating the MaxEntScan for the resulting sequence from the reference genome with the variant at the first position. This process is repeated as a sliding window across 9 5' bases and 23 3' splice site bases to obtain 2 arrays of MaxEntScan scores. The maximal value of each of these represents the window that contains the variant that is most likely to create a new splice site. For this maximal window, Introme calculates the corresponding MaxEntScan value using just the reference genome for comparison.

Introme is not aware of the directionality and strand of the gene within which each variant is located. As such, it performs the window approach described 4 times: forward, reverse, complement and reverse complement. The maximal value from all runs is taken as the most impactful.

MaxEntScan provides no way of interpreting the resulting score by default. Introme attempts to use the score to predict whether any given variant is likely to create a new splice site that overpowers the existing splice site. It does this by creating a "MaxEntScan Consequence" heuristic with NONE, LOW, MED and HIGH values for new splice site potential.

Score thresholds used to denote NONE, LOW, MED and HIGH consequences are at this point arbitrary but based on published work. First, https://doi.org/10.1002/humu.10295 reported that "The ideal MaxENT score is 11.81 for a 5' splice site and 13.59 for a 3' splice site. Theoretically, the larger the MaxENT value, the more efficient the splicing." Second, https://doi.org/10.1371/journal.pgen.1003613 contains distributions of MaxEnt scores across all introns in the genome in Figure 1. These show distributions centred on scores of 9 with a skew towards the right tail. From these, it appears that: a score of below 0 is almost never going to lead to a real splice site, a score of 0-4 is unlikely to lead to a real splice site, a score of 4-10 might lead to a new splice site and a score of above 10 is likely to lead to a new splice site.

Thus, these are the criteria for each of the consequence categories (they are applied to the 3' and 5' scores independently and the highest consequence is used in the annotation column):

* NONE: the window with the highest MaxEntScan score has a score lower than the corresponding reference genome OR the maximal value for the window with the variant is below 0
* LOW: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is 0-4 AND the reference genome score is >=4 smaller
* MED: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is 4-10 AND the reference genome score is >=4 smaller
* HIGH: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is >10 AND the reference genome score is >=6 smaller

The idea behind these heuristics is not to definitively indicate that a new splice site has been created, but instead to act as a guide. A more sophisticated analysis would be aware of the 5' and 3' scores for the real splice sites for the intron the variant is in and use this information intelligently with the sliding window approach and distance from splice site to predict impact. Of course, such an analysis should also be transcript and other upstream/downstream variant aware which adds another dimension of complexity...

