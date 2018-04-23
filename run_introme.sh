#!/bin/bash
cd /data/intronome/outputs
Fam=F1
ref_dir=~/data/intronome/references
intron=UCSC_intron.bed
5UTR=UCSC_5primeUTR.bed
3UTR=UCSC_3primeUTR.bed
2000US=references/UCSC_2000_upstream.bed
MotifFeatInput=homo_sapiens.GRCh37.motiffeatures.20161117.gff.gz
MotifFeat=Ens_MotifFeat
RegFeatInput=homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz 
RegFeat=Ens_RefFeat
input=~/data/intronome/data/R_170503_GINRAV_DNA_M001.hc.vqsr.vep.vcf.gz
daily_date=$(date +%d-%b)
echo 'beginning intersection' && \
bedtools intersect -header -a $input -b $ref_dir/$intron | bgzip > $Fam.$intron.vcf.gz && tabix -p vcf $Fam.$intron.vcf.gz && echo $intron 'done' \
&& bedtools intersect -header -a $input -b $ref_dir/5UTR | bgzip > $Fam.5UTR.vcf.gz && tabix -p vcf $Fam.5UTR.vcf.gz && echo $5UTR 'done' \
&& bedtools intersect -header -a $input -b $ref_dir/3UTR | bgzip > $Fam.3UTR.vcf.gz && tabix -p vcf out_dir/$Fam.5UTR.vcf.gz && echo $3UTR 'done' \
&& bedtools intersect -header -a $input -b $ref_dir/2000US | bgzip > $Fam.$2000US.vcf.gz && tabix -p vcf $Fam.$2000US.vcf.gz && echo $2000US 'done' \
&& bedtools intersect -header -a $input -b $ref_dir/MotifFeatInput | bgzip > $Fam.$MotifFeat.vcf.gz && tabix -p vcf $Fam.$MotifFeat.vcf.gz && echo $MotifFeatInput 'done' \
&& bedtools intersect -header -a $input -b $ref_dir/RegFeatInput | bgzip > $Fam.$RegFeat.vcf.gz && tabix -p vcf $Fam.$RegFeat.vcf.gz && echo $RegFeatInput 'done' \
&& echo 'all complete'

#Step 2- VCF anno. 

##IMPORTANT VCF FORMATTING NOTE####

##VCVFanno It requires bgzipped, tabix indexed BED, GFF VCF files as input. If your input VCF has incorectly formatted headers, VCFanno will not work. 
#eg ##INFO=<ID=CC,Number=1,Type=String,Description=""> This does not have any info between the "" and is invalid for VCFanno. Fix by removing headers with sed. E.g. 

#sed '1,14d' input.file > output.file && echo -e 'newheaderinfo' | cat - output.file > reheadered.output.file.

for VCF in $(ls *.vcf.gz)
do 
vcfanno conf.toml $VCF | bgzip > annotated.$VCF && tabix -p vcf annotated.$VCF && echo $VCF 'annotation complete'
done
#Use conf.toml file to specify what you want VCFanno to do. this is the one I used. 
#Step 3- Hard filtering your subsetted, annotated VCF
#Filtering thresholds will depend on your application. this is a guide. A useful doc is http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
for annotated_file in $(ls annotated.*.vcf.gz)
do 
zcat $annotated_file | perl -ne 'if ($_ =~ /gn_af=(.*?);/) { if ($1 < 0.5) { print $_; }}' | \
perl -ne 'if ($_ =~ /gn_ac=(.*?);/) { if ($1 <= 100) { print $_; }}' | \
perl -ne 'if ($_ =~ /DANN_score=(.*?);/) { if ($1 >= 0.8) { print $_; }}' | \
perl -ne 'if ($_ =~ /mgrb_ac=(.*?);/) { if ($1 <= 50) { print $_; }}' | \
perl -ne 'if ($_ =~ /fathmm-MKL_non-coding=(.*?);/) { if ($1 >= 0.6) { print $_; }}' | \
perl -ne 'if ($_ =~ /CADD_phred=(.*?);/) { if ($1 >= 10) { print $_; }}' | \
bgzip > filtered.$annotated_file && tabix filtered.$annotated_file
done


#step4- SHOW ME WHAT YOU GOT (ahem rick and morty). how many variants? ##includes headers!!##

wc -l filtered.*.vcf.gz > variant_count.$daily_date.txt
