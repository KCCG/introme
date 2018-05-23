#!/bin/bash

# Introme Bash Script
#
# The purpose of this script is to take a VCF file and first remove all variants that are
# not in regions of interest (intronic regions). Then, remaining variants are annotated
# with pathogenicity scores and global allele frequencies. Finally, hard filters are
# applied to annotations to remove variants that are almost definitely not pathogenic 
# (e.g. common in human populations).
#
# Developers:
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)
# Velimir Gayevskiy (vel@vel.nz)
#
# Usage: ./introme.sh input.vcf.gz family_ID inheritance_pattern
#
# TODO:
# * Need to let through variants that don't have an annotation rather than excluding them
#
# IMPORTANT VCF FORMATTING NOTE
# VCVFanno requires bgzipped, tabix indexed BED, GFF VCF files as input. If your input 
# VCF has incorectly formatted headers, VCFanno will not work. 
# e.g. ##INFO=<ID=CC,Number=1,Type=String,Description=""> This does not have any info 
# between the "" and is invalid for VCFanno. Fix by removing headers with sed. E.g. 
# sed '1,14d' input.file > output.file && echo -e 'newheaderinfo' | cat - output.file 
# > reheadered.output.file.

##################################

# Directories
working_dir=/data/intronome
ref_dir=/data/intronome/references
data_dir=/data/introme/data
out_dir/data/intronome/outputs

# Genome subsetting files (i.e. regions of interest)
intron=UCSC_intron.bed
5UTR=UCSC_5primeUTR.bed
3UTR=UCSC_3primeUTR.bed
2000US=references/UCSC_2000_upstream.bed
MotifFeatInput=homo_sapiens.GRCh37.motiffeatures.20161117.gff.gz
RegFeatInput=homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz 

# Script arguments
# Input VCF file
input_VCF=$1

# Family structure
Family=$2
inheritance_pattern=$3
PROBAND=D16-1922
MUM=D16-0862
DAD=D03-114
GRANNY=xxxxx
affected=$(echo %$PROBAND %$MUM %$DAD)
unaffected=$(echo %GRANNY)

# Hard filter cutoffs for each annotation
# Filtering thresholds will depend on your application, we urge you to carefully consider these 
gnomad_AF='<0.5'
gnomad_AC='<=100'
DANN_score='>=0.8'
MGRB_AC='<=50'
fathmm-MKL='>=0.6'
CADD_phred='>=10'

##################################
# STEP ONE: subsetting the VCF to genomic regions of interest

echo $(date +%x_%r) 'Beginning intersection'

bedtools intersect -header -a $input_VCF -b $ref_dir/$intron $ref_dir/$5UTR $ref_dir/$3UTR $ref_dir/$2000US $ref_dir/$MotifFeatInput $ref_dir/$RegFeatInput | bgzip > $out_dir/$Family.subset.vcf.gz
tabix -p vcf $out_dir/$Fam.subset.vcf.gz

echo $(date +%x_%r) 'Intersection complete'

##################################
# STEP TWO: annotate the subsetted VCF with useful information, to be used for filtering downstream
# The conf.toml file specifies what VCFanno should annotate

echo $(date +%x_%r) 'Beginning annotation'

vcfanno $ref_dir/conf.toml $out_dir/$Fam.subset.vcf.gz | \
bgzip > $out_dir/$Fam.subset.annotated.vcf.gz
tabix -p vcf $out_dir/$Fam.subset.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'
zcat $out_dir/annotated.$Fam.subset.vcf.gz | grep '#' | gzip > $out_dir/$Fam.subset.annotated.vcf.header.gz

##################################
# STEP THREE: Hard filtering the subsetted, annotated VCF

echo $(date +%x_%r) 'Beginning filtering' $out_dir/$Fam.subset.annotated.vcf.gz

bcftools filter -i'gn_af$gnomad_AF || gn_ac$gnomad_AC || DANN_score$DANN_score || mgrb_ac$MGRB_AC || fathmm-MKL_non-coding$fathmm-MKL || CADD_phred$CADD_phred' $out_dir/$Fam.subset.annotated.vcf.gz | bgzip > $out_dir/$Fam.subset.tmp.vcf.gz

# Add header information back again
zcat $out_dir/$Fam.subset.annotated.vcf.header.gz zcat $out_dir/$Fam.subset.tmp.vcf.gz | \
bgzip > $out_dir/$Fam.subset.annotated.filtered.vcf.gz
tabix $out_dir/$Fam.subset.annotated.filtered.vcf.gz
rm $out_dir/$Fam.subset.tmp.vcf.gz

echo $(date +%x_%r) 'Filtering complete' 

##################################
# STEP 4: family filtering

if [ "$inheritance_pattern" = "denovo"]; then
	echo "applying de novo genotype filter"
	bcftools filter -i'0/1' -f'%$PROBAND\n' || -e'0/1' -f'%$MUM %$DAD\n' || -e'1/1' -f'%$MUM %DAD/n' $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz
elif ["$inheritance_pattern" = "recessive"]; then
	echo "applying homozygous recessive genotype filter"
	bcftools filter -i'1/1' -f'%$PROBAND\n' || -i'0/1' -f '%$MUM %$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz
elif [ "$inheritance_pattern" = "dominant"]; then
	echo "applying dominant genotype filter"
	bcftools filter -i'0/1' -f'$affected\n' || -e'1/1' -f'$affected\n' || -e'1/1' -f'$unaffected\n' || -e'0/1' -f'$unaffected\n'  $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz
# elif [ "$inheritance_pattern" = "comphet"]; then
	# echo "applying compound heterozygous genotype filter"
	# bcftools filter -i'0/1' -f'%$PROBAND %$MUM\n' || -e'1/1' -f'%$PROBAND %$MUM\n' || -e'0/1' -f'%$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$MUM.variants.txt
	# bcftools filter -i'0/1' -f'%$PROBAND %$DAD\n' || -e'1/1' -f'%$PROBAND %$DAD\n' || -e'0/1' -f'%$MUM\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$DAD.variants.txt
fi

##################################

exit