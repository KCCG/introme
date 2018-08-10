# What is Introme?

Introme identifies genetic variants that are likely to disrupt existing splice sites or create new ones. Currently, these variants are discarded as 'low impact' in research and diagnostic settings, but have the potential to be just as damaging as exonic variants in disrupting normal gene expression.

Introme accepts a joint-called VCF file containing up to millions of variants found in a typical family trio. It annotates the variants with a number of prediction scores, allele frequencies and pathogenicity estimates, then applies a number of filters to enrich for variants with pathogenic potential.

# Developers

Introme was conceived and initially specced as an idea by Dr. Mark Cowley and Dr. Velimir Gayevskiy. Sarah Beecroft carried out initial implementation. Dr. Velimir Gayevskiy reimplemented and added many features.

# Requirements

## You need the following software installed

Introme was developed on macOS and runs with the following dependencies which can be mostly found using the Homebrew package manager for macOS. Since it's a Bash script, it should work just fine on Linux with the same dependencies installed.

* bedtools
* tabix
* vcfanno
* bcftools
* bgzip
* MaxEntScan (supplied with Introme from http://genes.mit.edu/burgelab/maxent/download/)

## The VCF file supplied should be

* Created using GATK HaplotypeCaller (other variant callers may work but are currently untested)
* Run through GATK VQSR (again, not strictly required)
* Decomposed and normalised (i.e. no multi-allelic variants) with vt
* Annotated with VEP (optional but very useful for knowing the standard predictions for each variant with transcripts)

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

* `-a` -- one or more affected sample names (specify multiple times if more than one), must have the same sample name as in the VCF file
* `-b` -- path to search space BED file
* `-i` -- inheritance pattern, can be one of "heterozygous", "homozygous" or "none"
* `-p` -- output prefix
* `-v` -- path to VCF file
* `-r` -- path to reference genome fasta file (must be indexed with samtools)

## Example

./run\_introme.sh -b subsetting/gencode.v28lift37.annotation.gtf.bed.gz -v input/Fam1_jointcall.hc.vqsr.decomposed.normalised.vep.vcf.gz -p Fam1 -i heterozygous -a A001C -a A002B -r annotations/hs37d5.fasta-index/genome.fa

# MaxEntScan consequence logic

MaxEntScan produces a numeric score for a given string of 9 bases for 5' splice sites and 23 bases for 3' splice sites. This score denotes the affinity of the bases for acting as a splice site. However, in the context of determining whether a single SNV in an intron is creating a new splice site, it is not obvious whether the splice site created is a 5' or 3' splice site and also where within the splice site the variant falls.

Introme attempts to find the maximal splice site creation potential of the variant. It does this by first supposing that the variant is the first base of a new splice site and calculating the MaxEntScan for the resulting sequence from the reference genome with the variant at the first position. This process is repeated as a sliding window across 9 5' bases and 23 3' splice site bases to obtain 2 arrays of MaxEntScan scores. The maximal value of each of these represents the window that contains the variant that is most likely to create a new splice site. For this maximal window, Introme calculates the corresponding MaxEntScan value using just the reference genome for comparison.

Genes can be located on the positive (forward) or negative (reverse) strands. If a gene is on a positive (forward) strand, it will be read in the same way as the reference genome sequence. If it's on the negative (reverse) strand, it will be read in reverse complement. Introme is aware of the direction of the gene(s) within which the variant falls. If there is more than one and all genes are pointing in the same direction, MaxEntScan will be calculated only in that direction (determining the reverse complement when necessary). Otherwise, it will be calculated in both directions and the maximum value will be reported along with the direction it was found in.

MaxEntScan provides no way of interpreting the resulting score by default. Introme attempts to use the score to predict whether any given variant is likely to create a new splice site that overpowers the existing splice site. It does this by creating a "MaxEntScan Consequence" heuristic with NONE, LOW, MED and HIGH values for new splice site potential.

Score thresholds used to denote NONE, LOW, MED and HIGH consequences are at this point arbitrary but based on published work. First, https://doi.org/10.1002/humu.10295 reported that "The ideal MaxENT score is 11.81 for a 5' splice site and 13.59 for a 3' splice site. Theoretically, the larger the MaxENT value, the more efficient the splicing." Second, https://doi.org/10.1371/journal.pgen.1003613 contains distributions of MaxEnt scores across all introns in the genome in Figure 1. These show distributions centred on scores of 9 with a skew towards the right tail. From these, it appears that: a score of below 0 is almost never going to lead to a real splice site, a score of 0-4 is unlikely to lead to a real splice site, a score of 4-10 might lead to a new splice site and a score of above 10 is likely to lead to a new splice site.

Thus, these are the criteria for each of the consequence categories (they are applied to the 3' and 5' scores independently and the highest consequence is used in the annotation column):

* NONE: the window with the highest MaxEntScan score has a score lower than the corresponding reference genome OR the maximal value for the window with the variant is below 0
* LOW: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is 0-4 AND the reference genome score is >=4 smaller
* MED: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is 4-10 AND the reference genome score is >=4 smaller
* HIGH: the window with the highest MaxEntScan score has a score higher than the corresponding reference genome AND the score is >10 AND the reference genome score is >=6 smaller

The idea behind these heuristics is not to definitively indicate that a new splice site has been created, but instead to act as a guide. A more sophisticated analysis would be aware of the 5' and 3' scores for the real splice sites for the intron the variant is in and use this information intelligently with the sliding window approach and distance from splice site to predict impact. Of course, such an analysis should also be transcript and other upstream/downstream variant aware which adds another dimension of complexity...

# Prediction of damage to existing splice sites and branch points

The SPIDEX and dbscSNV annotations are used to predict whether a variant is damaging existing splice sites. SPIDEX looks at all exonic variants and 300bp into the introns while dbscSNV looks in the splicing regions only.

The Branchpointer annotation is used to finding variants that are disrupting branch points and potentially leading to a failure of splicing. You will need to manually determine whether there are other branch points nearby that may take over for the damaged one. 

# Interpretation of results

You will need to subset the results file from Introme for each annotation: SPIDEX, dbscSNV, Branchpointer and MaxEntScan. The underlying genetic mechanisms are distinct so you will need to subset multiple times to interrogate each way splicing can be disrupted (by damage to an existing splice site, by damage to a branch point and by the creating of a new splice site). Any interesting variants will need to be interpreted in light of the sequence around them. For example, finding a variant in an intron that is creating a new 5' splice site means we are hypothesising that this variant is at the very end of a new exon, therefore we need to look at up to 1000bp upstream for an existing 3' splice site and 18-44bp upstream of that for a branch point.

Do that like so:

1) Pull out the reference genome sequence up/downstream of the variant depending on if it's a 3' or a 5' splice site-creating variant.

samtools faidx ./genome.fa 4:27993059-27993759 

Reverse complement this if needed (if the gene is on the reverse strand)

echo "<sequence>" | tr "[ATCG]" "[TAGC]" | rev

2) Walk through the sequence as windows and run MaxEntScan on each window.

for i in {0..588}; do echo ${seq:(( 0 + $i )):9}; done | perl ./MaxEntScan/score5.pl - > ./5_prime_output.txt

Find largest values (in Excel) and find the corresponding windows in the sequence pulled out in step 1.

3) Functionally validate the finding. Good luck with that one.