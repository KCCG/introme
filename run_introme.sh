#!/bin/bash

# Introme Bash Script
#
# Developers:
# Velimir Gayevskiy (vel@vel.nz)
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)




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
CADD_phred='>=10'
min_QUAL='>=200' # The QUAL VCF field
max_DP='>=10' # The sample with the highest depth in the VCF must have a higher depth than this value

##################################
# STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)

echo $(date +%x_%r) 'Beginning subsetting'
echo $(gzip -d -c $input_VCF | wc -l) 'VCF lines prior to subsetting'

bedtools intersect -header -u -a $input_VCF -b $input_BED | bgzip > $out_dir/$prefix.subset.vcf.gz # -u for unique record in VCF, otherwise multiple variants are output for overlapping introns
tabix -p vcf $out_dir/$prefix.subset.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.vcf.gz | wc -l) 'VCF lines after subsetting'
echo $(date +%x_%r) 'Subsetting complete'

##################################
# STEP 2: familial filtering (second because it gets rid of the second most variants)

echo $(date +%x_%r) 'Beginning familial filtering'

#if [ "$inheritance_pattern" = "denovo"]; then
#	echo "applying de novo genotype filter"
	bcftools filter --threads 8 -i"FORMAT/GT[0]='0/0' && FORMAT/GT[1]='0/0' && FORMAT/GT[2]='0/1'" $out_dir/$prefix.subset.vcf.gz | bgzip > $out_dir/$prefix.subset.denovo.vcf.gz
	tabix $out_dir/$prefix.subset.denovo.vcf.gz
	echo $(gzip -d -c $out_dir/$prefix.subset.denovo.vcf.gz | wc -l) 'VCF lines after de novo filter'
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

echo $(date +%x_%r) 'Familial filtering complete'

##################################
# STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream

echo $(date +%x_%r) 'Beginning annotation'

vcfanno -p 8 conf.toml $out_dir/$prefix.subset.denovo.vcf.gz | bgzip > $out_dir/$prefix.subset.denovo.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
tabix -p vcf $out_dir/$prefix.subset.denovo.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'

##################################
# STEP 4: Hard filtering the subsetted, annotated VCF

echo $(date +%x_%r) 'Beginning filtering'

bcftools filter --threads 8 -i"FILTER='PASS' && TYPE='snp' && QUAL$min_QUAL && MAX(DP)$max_DP && (mgrb_af$MGRB_AF || mgrb_af='.') && (gn_pm_af$gnomad_AF || gn_pm_af='.') && (cadd_phred$CADD_phred || cadd_phred='.')" $out_dir/$prefix.subset.denovo.annotated.vcf.gz | bgzip > $out_dir/$prefix.subset.denovo.annotated.filtered.vcf.gz
tabix -p vcf $out_dir/$prefix.subset.denovo.annotated.filtered.vcf.gz

#bcftools filter -i'gn_af$gnomad_AF || gn_ac$gnomad_AC || DANN_score$DANN_score || mgrb_ac$MGRB_AC || fathmm-MKL_non-coding$fathmm-MKL || CADD_phred$CADD_phred' $out_dir/$Fam.subset.annotated.vcf.gz | bgzip > $out_dir/$Fam.subset.tmp.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.denovo.annotated.filtered.vcf.gz | wc -l) 'VCF lines after filtering'

echo $(date +%x_%r) 'Filtering complete' 

##################################

# For MaxEntScan:
# echo "ctctactactatctatctagatc" | perl score3.pl - | cut -f 2

##################################

exit