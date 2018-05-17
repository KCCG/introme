#!/bin/bash
#usage is ./introme.sh input.joint.vcf family_ID inheritance_pattern
#input joint VCF
input_VCF=$1
#used R_170503_GINRAV_DNA_M001.hc.vqsr.vep.vcf.gz
#Family structure
Family=$2
inheritance_pattern=$3
PROBAND=D16-1922
MUM=D16-0862
DAD=D03-114
GRANNY=xxxxx
affected=$(echo %$PROBAND %$MUM %$DAD)
unaffected=$(echo %GRANNY)
#directories
working_dir=/data/intronome
out_dir/data/intronome/outputs
ref_dir=/data/intronome/references
data_dir=/data/introme/data
#genome subset files
intron=UCSC_intron.bed
5UTR=UCSC_5primeUTR.bed
3UTR=UCSC_3primeUTR.bed
2000US=references/UCSC_2000_upstream.bed
MotifFeatInput=homo_sapiens.GRCh37.motiffeatures.20161117.gff.gz
RegFeatInput=homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz 
#Hard filter cutoffs for step 3
gnomad_AF='<0.5'
gnomad_AC='<=100'
DANN_score='>=0.8'
MGRB_AC='<=50'
fathmm-MKL='>=0.6'
CADD_phred='>=10

###STEP ONE- subsetting your joint VCF to genomic regions of interest
echo $(date +%x_%r) 'beginning intersection'
bedtools intersect -header -a $input_VCF -b $ref_dir/$intron $ref_dir/$5UTR $ref_dir/$3UTR $ref_dir/$2000US $ref_dir/$MotifFeatInput $ref_dir/$RegFeatInput | bgzip > $out_dir/$Family.subet.vcf.gz && tabix -p vcf $out_dir/$Fam.subset.vcf.gz
echo $(date +%x_%r) 'intersection complete'

###STEP TWO- annotating your subsetted VCF/s with useful information, to be used for filtering downstream
#Use conf.toml file to specify what you want VCFanno to do (example at https://github.com/SarahBeecroft/introme)

                                                      ##IMPORTANT VCF FORMATTING NOTE####
##VCVFanno It requires bgzipped, tabix indexed BED, GFF VCF files as input. If your input VCF has incorectly formatted headers, VCFanno will not work. 
#eg ##INFO=<ID=CC,Number=1,Type=String,Description=""> This does not have any info between the "" and is invalid for VCFanno. Fix by removing headers with sed. E.g. 
#sed '1,14d' input.file > output.file && echo -e 'newheaderinfo' | cat - output.file > reheadered.output.file.

echo $(date +%x_%r) 'annotation beginning'
vcfanno $ref_dir/conf.toml $out_dir/$Fam.subset.vcf.gz | \
bgzip > $out_dir/$Fam.subset.annotated.vcf.gz
tabix -p vcf $out_dir/$Fam.subset.annotated.vcf.gz
echo $(date +%x_%r) 'annotation complete'
zcat $out_dir/annotated.$Fam.subset.vcf.gz | grep '#' | gzip > $out_dir/$Fam.subset.annotated.vcf.header.gz

#Step 3- Hard filtering your subsetted, annotated VCF
#Filtering thresholds will depend on your application. this is a guide. A useful doc is http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
echo $(date +%x_%r) 'filtering' $out_dir/$Fam.subset.annotated.vcf.gz
bcftools filter -i'gn_af$gnomad_AF || gn_ac$gnomad_AC || DANN_score$DANN_score || mgrb_ac$MGRB_AC || fathmm-MKL_non-coding$fathmm-MKL || CADD_phred$CADD_phred' $out_dir/$Fam.subset.annotated.vcf.gz | bgzip > $out_dir/$Fam.subset.tmp.vcf.gz

#add header information back again
zcat $out_dir/$Fam.subset.annotated.vcf.header.gz zcat $out_dir/$Fam.subset.tmp.vcf.gz | \
bgzip > $out_dir/$Fam.subset.annotated.filtered.vcf.gz
tabix $out_dir/$Fam.subset.annotated.filtered.vcf.gz
rm $out_dir/$Fam.subset.tmp.vcf.gz
echo $(date +%x_%r) 'completed filtering' 

##STEP 4 - family filtering
#de novo filter 
if [ "$inheritance_pattern" = "de novo"]; then
  echo "applying de novo genotype filter"
  bcftools filter -i'0/1' -f'%$PROBAND\n' || -e'0/1' -f'%$MUM %$DAD\n' || -e'1/1' -f'%$MUM %DAD/n' $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.de_novo.vcf.gz
  elif ["$inheritance_pattern" = "homozygous recessive"]; then
  echo "applying homozygous recessive genotype filter"
  bcftools filter -i'1/1' -f'%$PROBAND\n' || -i'0/1' -f '%$MUM %$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.hom_rec.vcf.gz
  elif [ "$inheritance_pattern" = "dominant"]; then
  echo "applying dominant genotype filter"
bcftools filter -i'0/1' -f'$affected\n' || -e'1/1' -f'$affected\n' || -e'1/1' -f'$unaffected\n' || -e'0/1' -f'$unaffected\n'  $out_dir/$Fam.subset.annotated.filtered.vcf.gz | bgzip > $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz && tabix $out_dir/$Fam.subset.annotated.filtered.dominant.vcf.gz
 # elif [ "$inheritance_pattern" = "compound heterozygous"]; then
  # echo "applying compound heterozygous genotype filter"
  # bcftools filter -i'0/1' -f'%$PROBAND %$MUM\n' || -e'1/1' -f'%$PROBAND %$MUM\n' || -e'0/1' -f'%$DAD\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$MUM.variants.txt
  # bcftools filter -i'0/1' -f'%$PROBAND %$DAD\n' || -e'1/1' -f'%$PROBAND %$DAD\n' || -e'0/1' -f'%$MUM\n' out_dir/$Fam.subset.annotated.filtered.vcf.gz > $out_dir/$DAD.variants.txt
fi
