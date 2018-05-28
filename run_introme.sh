#!/bin/bash

# Introme Bash Script
#
# Developers:
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)
# Velimir Gayevskiy (vel@vel.nz)
#
# TODO:
# * Need to let through variants that don't have an annotation rather than excluding them
#



##################################
# DEFINE VARIABLES

# Directories
anno_dir=./annotations
out_dir=./output

# Script arguments
# Input VCF file
input_VCF=$1

# Input BED file (i.e. regions of interest)
input_BED=$2

# Output file prefix
prefix=$3

# Family structure
#inheritance_pattern=$4
#PROBAND=D16-1922
#MUM=D16-0862
#DAD=D03-114
#GRANNY=xxxxx
#affected=$(echo %$PROBAND %$MUM %$DAD)
#unaffected=$(echo %GRANNY)

# Hard filter cutoffs for each annotation
# Filtering thresholds will depend on your application, we urge you to carefully consider these 
MGRB_AF='<=0.01'
gnomad_AF='<=0.01'

#DANN_score='>=0.8'
#fathmm_MKL='>=0.6'
#CADD_phred='>=10'

##################################
# STEP ONE: subsetting the VCF to genomic regions of interest

echo $(date +%x_%r) 'Beginning subsetting'
echo $(gzip -d -c $input_VCF | wc -l) 'VCF lines prior to subsetting'

bedtools intersect -header -u -a $input_VCF -b $input_BED | bgzip > $out_dir/$prefix.subset.vcf.gz # -u for unique record in VCF, otherwise multiple variants are output for overlapping introns
tabix -p vcf $out_dir/$prefix.subset.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.vcf.gz | wc -l) 'VCF lines after subsetting'
echo $(date +%x_%r) 'Subsetting complete'

##################################
# STEP TWO: annotate the subsetted VCF with useful information, to be used for filtering downstream

echo $(date +%x_%r) 'Beginning annotation'

vcfanno -p 4 conf.toml $out_dir/$prefix.subset.vcf.gz | bgzip > $out_dir/$prefix.subset.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
tabix -p vcf $out_dir/$prefix.subset.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'

#zcat $out_dir/annotated.$Fam.subset.vcf.gz | grep '#' | gzip > $out_dir/$Fam.subset.annotated.vcf.header.gz

##################################
# STEP THREE: Hard filtering the subsetted, annotated VCF

echo $(date +%x_%r) 'Beginning filtering'

bcftools filter -i"mgrb_af$MGRB_AF && gn_af$gnomad_AF" $out_dir/$prefix.subset.annotated.vcf.gz | bgzip > $out_dir/$prefix.subset.annotated.filtered.vcf.gz

# TODO which gnomad column is best to use?

#bcftools filter -i'gn_af$gnomad_AF || gn_ac$gnomad_AC || DANN_score$DANN_score || mgrb_ac$MGRB_AC || fathmm-MKL_non-coding$fathmm-MKL || CADD_phred$CADD_phred' $out_dir/$Fam.subset.annotated.vcf.gz | bgzip > $out_dir/$Fam.subset.tmp.vcf.gz

# Add header information back again
#zcat $out_dir/$Fam.subset.annotated.vcf.header.gz zcat $out_dir/$Fam.subset.tmp.vcf.gz | \
#bgzip > $out_dir/$Fam.subset.annotated.filtered.vcf.gz
#tabix $out_dir/$Fam.subset.annotated.filtered.vcf.gz
#rm $out_dir/$Fam.subset.tmp.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.annotated.filtered.vcf.gz | wc -l) 'VCF lines after filtering'

echo $(date +%x_%r) 'Filtering complete' 

##################################
# STEP 4: family filtering

# TODO apply de novo inheritance pattern and see how many are left!
# Read man page for bcftools for how best to do this

#if [ "$inheritance_pattern" = "denovo"]; then
#	echo "applying de novo genotype filter"
#	bcftools filter -i'0/1' -f'%$PROBAND\n' || -e'0/1' -f'%$MUM %$DAD\n' || -e'1/1' -f'%$MUM %DAD/n' $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz
#elif ["$inheritance_pattern" = "recessive"]; then
#	echo "applying homozygous recessive genotype filter"
#	bcftools filter -i'1/1' -f'%$PROBAND\n' || -i'0/1' -f '%$MUM %$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz
#elif [ "$inheritance_pattern" = "dominant"]; then
#	echo "applying dominant genotype filter"
#	bcftools filter -i'0/1' -f'$affected\n' || -e'1/1' -f'$affected\n' || -e'1/1' -f'$unaffected\n' || -e'0/1' -f'$unaffected\n'  $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz
# elif [ "$inheritance_pattern" = "comphet"]; then
	# echo "applying compound heterozygous genotype filter"
	# bcftools filter -i'0/1' -f'%$PROBAND %$MUM\n' || -e'1/1' -f'%$PROBAND %$MUM\n' || -e'0/1' -f'%$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$MUM.variants.txt
	# bcftools filter -i'0/1' -f'%$PROBAND %$DAD\n' || -e'1/1' -f'%$PROBAND %$DAD\n' || -e'0/1' -f'%$MUM\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$DAD.variants.txt
#fi

# For MaxEntScan:
# echo "ctctactactatctatctagatc" | perl score3.pl - | cut -f 2

##################################

exit