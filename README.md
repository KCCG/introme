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
* MaxEntScan (supplied from http://genes.mit.edu/burgelab/maxent/download/)

## The VCF file supplied should be

* Created using GATK HaplotypeCaller (other variant callers may work but are currently untested)
* Decomposed and normalised (i.e. no multi-allelic variants) with vt
* Annotated with VEP (optional but very useful for knowing what genes the variants are in)

## Annotations

* CADD v1.3 VCF created using the instructions at: https://github.com/brentp/vcfanno/blob/master/docs/examples/cadd.md
* gnomad.genomes.sites.merged.AF\_AC_AN\_only.vcf.gz
* MGRB variant allele frequencies (these are not public so remove this annotation manually)

# Using Introme

## Options

* -a -- one or more affected sample names (specify multiple times if more than one), must have the same sample name as in the VCF file
* -b -- path to search space BED file
* -i -- inheritance pattern, can be one of "denovo" or "autrec"
* -p -- output prefix
* -v -- path to VCF file
* -r -- path to reference genome fasta file (must be indexed)

## Example

./run\_introme.sh -b subsetting/UCSC\_introns.bed.gz -v input/Fam1_jointcall.hc.vqsr.decomposed.normalised.vep.vcf.gz -p Fam1 -i denovo -a A001C -r annotations/hs37d5.fasta-index/genome.fa
