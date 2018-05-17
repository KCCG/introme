#!/bin/bash
#input joint VCF
input_VCF=$1
#used R_170503_GINRAV_DNA_M001.hc.vqsr.vep.vcf.gz
#Family structure
#define affected samples
PROBAND=D16-1922
AFF2=D16-0862
AFF3=D03-114
#define UNaffected samples
#HEALTHY1=
#HEALTHY2=
#directories
working_dir=/data/intronome/outputs
ref_dir=/data/intronome/references
data_dir=/data/introme/data
#genome subset files
intron=UCSC_intron.bed
5UTR=UCSC_5primeUTR.bed
3UTR=UCSC_3primeUTR.bed
2000US=references/UCSC_2000_upstream.bed
MotifFeatInput=homo_sapiens.GRCh37.motiffeatures.20161117.gff.gz
MotifFeat=Ens_MotifFeat
RegFeatInput=homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz 
RegFeat=Ens_RefFeat
#Patient selection
Family=F1
#Hard filter cutoffs for step 3
gnomad_AF= <0.5
gnomad_AC= <=100
DANN_score= >=0.8
MGRB_AC= <=50
fathmm-MKL= >=0.6
CADD_phred= >=10

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
&& tabix -p vcf $out_dir/$Fam.subset.annotated.vcf.gz
echo $(date +%x_%r) 'annotation complete'
zcat $out_dir/annotated.$Fam.subset.vcf.gz | grep '#' | gzip > $out_dir/$Fam.subset.annotated.vcf.header.gz

#Step 3- Hard filtering your subsetted, annotated VCF
#Filtering thresholds will depend on your application. this is a guide. A useful doc is http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
echo $(date +%x_%r) 'filtering' $out_dir/$Fam.subset.annotated.vcf.gz

zcat $out_dir/$Fam.subset.annotated.vcf.gz | perl -ne 'if ($_ =~ /gn_af=(.*?);/) { if ($1 $gnomad_AF) { print $_; }}' | \
perl -ne 'if ($_ =~ /gn_ac=(.*?);/) { if ($1 $gnomad_AC) { print $_; }}' | \
perl -ne 'if ($_ =~ /DANN_score=(.*?);/) { if ($1 $DANN_score) { print $_; }}' | \
perl -ne 'if ($_ =~ /mgrb_ac=(.*?);/) { if ($1 $MGRB_AC) { print $_; }}' | \
perl -ne 'if ($_ =~ /fathmm-MKL_non-coding=(.*?);/) { if ($1 $fathmm-MKL) { print $_; }}' | \
perl -ne 'if ($_ =~ /CADD_phred=(.*?);/) { if ($1 $CADD_phred) { print $_; }}' | \
bgzip > $out_dir/$Fam.subset.tmp.vcf.gz
#add header information back again
zcat $out_dir/$Fam.subset.annotated.vcf.header.gz zcat $out_dir/$Fam.subset.tmp.vcf.gz | \
bgzip > $out_dir/$Fam.subset.annotated.filteredvcf.gz
tabix $out_dir/$Fam.subset.annotated.filteredvcf.gz
echo $(date +%x_%r) 'completed filtering' $annotatedVCF

##STEP 4 - family filtering
#de novo filter 
PROBAND=D16-1922
AFF2=D16-0862
AFF3=D03-114

Dominant
bcftools filter -i'0/1' -f'%$PROBAND %$AFF2 %$AFF3\n' $file > dominant.gt_filter.%file

De novo
bcftools filter -i'0/1' -f'%$PROBAND\n' || -e'0/1' -f'%$MUM %$DAD\n' $file > denovo.gt_filter.%file

Homo recessive
bcftools filter -i'1/1' -f'%$PROBAND %$AFF2 %$AFF3\n' || -i'0/1' -f '%$MUM %$DAD\n' $file > homRec.gt_filter.%file
if filter=de_novo
then
do

print where GT 0/1
