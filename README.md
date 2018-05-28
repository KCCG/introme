# What is Introme?

The purpose of this script is to take a VCF file and first remove all variants that are not in regions of interest (intronic regions). Then, remaining variants are annotated with pathogenicity scores and global allele frequencies. Finally, hard filters are applied to annotations to remove variants that are almost definitely not pathogenic (e.g. common in human populations).

# Using Introme

./introme.sh input.vcf.gz UCSC\_introns.bed.gz output\_prefix

# Requirements

You need the following software installed:
* bedtools
* tabix
* vcfanno
* bcftools

The VCF file supplied should be:
* Decomposed and normalised (i.e. no multi-allelic variants)

# vcfanno annotations

vcfanno requires bgzipped, tabix indexed BED, GFF or VCF files as input. If your input VCF has incorectly formatted headers, vcfanno will not work. e.g. ##INFO=<ID=CC,Number=1,Type=String,Description=""> This does not have any info between the "" and is invalid for vcfanno. Fix by removing headers with sed. E.g. sed '1,14d' input.file > output.file && echo -e 'newheaderinfo' | cat - output.file > reheadered.output.file.