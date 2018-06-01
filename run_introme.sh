#!/bin/bash

# Introme Bash Script
#
# Developers:
# Velimir Gayevskiy (vel@vel.nz)
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)
#
# TODO:
# 1) Perhaps pre-calculate all possible max MaxEntScan scores/sequences for every possible SNP in the reference genome? Back-of-the-envelope estimate for this is 50,000 days on 1 CPU (for just 5' though) so with access to a few thousand it wouldn't take long.





##################################
# CUSTOM FUNCTIONS

# Check whether a string is an array element
# Usage: stringinarray "string" "${array[@]}"
# $? == 0 if true, otherwise 1
stringinarray () {
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}

##################################
# DEFINE VARIABLES

# Directories
anno_dir='./annotations'
out_dir='./output'

# Hard filter cutoffs for each annotation, filtering thresholds will depend on your application, we urge you to carefully consider these 
MGRB_AF='<=0.01' # Minimum MGRB allele frequency (healthy Australian population)
gnomad_popmax_AF='<=0.01' # Minimum gnomAD allele frequency in the population with the highest allele frequency in gnomAD
CADD_phred='>=10' # Minimum CADD score (phred-scaled)
min_QUAL='>=200' # The QUAL VCF field
max_DP='>=20' # The sample with the highest depth in the VCF must have a higher depth than this value

##################################
# INTROME ARGUMENTS

while getopts "a:b:v:p:i:r:" opt; do
    case $opt in
        a) affected_samples+=("$OPTARG");; # Array of affected samples
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        v) input_VCF="$OPTARG";; # Input VCF file
        p) prefix="$OPTARG";; # Output file prefix
        i) inheritance_pattern="$OPTARG";; # Inheritance pattern to use
        r) reference_genome="$OPTARG";; # Path to the reference genome used for mapping
    esac
done
shift $((OPTIND -1)) # Remove parsed options and args from $@ list

# Make sure each required argument has a non-empty value
if [ -z $affected_samples ]; then
	echo "No affected samples supplied."
	exit 1
elif [ -z $input_BED ]; then
	echo "No restriction BED file supplied."
	exit 1
elif [ -z $input_VCF ]; then
	echo "No path to an input VCF file has been supplied."
	exit 1
elif [ -z $prefix ]; then
	echo "No output prefix has been supplied."
	exit 1
elif [ -z $inheritance_pattern ]; then
	echo "No inheritance pattern has been supplied."
	exit 1
elif [ -z $reference_genome ]; then
	echo "No reference genome file has been supplied."
	exit 1
fi

##################################
# STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)

echo $(date +%x_%r) 'Beginning subsetting'
echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants prior to subsetting'

bedtools intersect -header -u -a $input_VCF -b $input_BED | bgzip > $out_dir/$prefix.subset.vcf.gz # -u for unique record in VCF, otherwise multiple variants are output for overlapping introns
tabix -p vcf $out_dir/$prefix.subset.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.vcf.gz | grep -v '^#' | wc -l) 'variants after subsetting'
echo $(date +%x_%r) 'Subsetting complete'

##################################
# STEP 2: familial filtering (second because it gets rid of the second most variants)

echo $(date +%x_%r) 'Beginning familial filtering'

# Capture each sample name in the VCF into an array
vcfsamples=( $(bcftools query -l $out_dir/$prefix.subset.vcf.gz) )

# Variable to hold the bcftools genotype filter string
genotype_filters=""

# Variable to hold the current VCF sample index, used in referring to the right sample in filtering
i=0

# Go through every VCF sample name (they are ordered as they are in the VCF which is important for the index)
while [[ $i -le $((${#vcfsamples[@]}-1)) ]]; do 
	# Start of the genotype filter string
	genotype_filters+=" && FORMAT/GT[$i]='"
		
		# Run custom function to check if the current VCF sample has been supplied as an affected sample
		stringinarray "${vcfsamples[$i]}" "${affected_samples[@]}"
		
		# If the current VCF sample is affected
		if [ $? == 0 ]; then
			if [[ $inheritance_pattern == "denovo" ]]; then
				genotype_filters+="0/1"
			elif [[ $inheritance_pattern == "autrec" ]]; then
				genotype_filters+="1/1"
			fi
		# If the current VCF sample is not affected
		else
			if [[ $inheritance_pattern == "denovo" ]]; then
				genotype_filters+="0/0"
			elif [[ $inheritance_pattern == "autrec" ]]; then
				genotype_filters+="0/1"
			fi
		fi

	# Close off genotype filter string
	genotype_filters+="'"

	# Iterate the $i variable
	((i++))
done

# Cut the 4 first characters off the filter string (the " && ")
genotype_filters=$(echo $genotype_filters | cut -c 4-)

bcftools filter --threads 8 -i"$genotype_filters" $out_dir/$prefix.subset.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.vcf.gz
tabix $out_dir/$prefix.subset.inheritancefilter.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.inheritancefilter.vcf.gz | grep -v '^#' | wc -l) 'variants after familial filter'
echo $(date +%x_%r) 'Familial filtering complete'

##################################
# STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream

echo $(date +%x_%r) 'Beginning annotation'

vcfanno -p 8 conf.toml $out_dir/$prefix.subset.inheritancefilter.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
tabix -p vcf $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'

##################################
# STEP 4: Hard filtering the subsetted, annotated VCF

echo $(date +%x_%r) 'Beginning filtering'

bcftools filter --threads 8 -i"FILTER='PASS' && TYPE='snp' && QUAL$min_QUAL && MAX(DP)$max_DP && (mgrb_af$MGRB_AF || mgrb_af='.') && (gn_pm_af$gnomad_popmax_AF || gn_pm_af='.') && (cadd_phred$CADD_phred || cadd_phred='.')" $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz
tabix -p vcf $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Filtering complete'

##################################
# STEP 5: Output final list of variants as a spreadsheet-friendly TSV

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/cadd_phred\t%INFO/mgrb_af\t%INFO/gn_pm_af\t%INFO/CSQ[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz > $out_dir/$prefix.introme.tsv

##################################
# STEP 6: Calculate the maximum MaxEntScan value for sliding windows around each variant

echo $(date +%x_%r) 'Beginning MaxEntScan calculation'

# Go through each line of the results TSV
cat $out_dir/$prefix.introme.tsv | \
while read line; do
	# If the line if the header line
	if [ "${line:0:1}" == "#" ]; then
		# Add to the header line
		echo "$line"$'\t'"5'_max_MaxEntScan_value"$'\t'"5'_max_MaxEntScan_sequence"$'\t'"5'_max_MaxEntScan_coordinates"$'\t'"5'_max_MaxEntScan_germline_value"$'\t'"5'_max_MaxEntScan_germline_sequence" > $out_dir/$prefix.introme.annotated.tsv
		
		# Don't do any further processing on this line
		continue
	fi
	
	max_maxentscan=''
	max_maxentscan_sequence=''
	max_maxentscan_coordinates=''
	current_chr=$(echo "$line" | cut -f 1)
	current_pos=$(echo "$line" | cut -f 2)
	current_alt=$(echo "$line" | cut -f 4)
	
	# $i is the offset for determining the coordinate range
	for i in {0..8}; do
		# Determine the putative splice site coordinate range
		current_coordinate_range="$current_chr:"$(($current_pos - 8 + $i))"-"$(($current_pos + i))
		
		# Extract the putative splice site sequence from the reference genome
		current_splice_site=$(samtools faidx $reference_genome $current_coordinate_range | grep -v '^>')
		
		# Mutate the correct base in the extracted splice site to match the current variant
		current_splice_site=${current_splice_site:0:(( 8 - $i ))}$current_alt${current_splice_site:(( 8 - $i + 1 ))}
		
		# Calculate MaxEntScan value
		current_maxentscan=$(echo $current_splice_site | perl MaxEntScan/score5.pl - | cut -f 2)
		
		# If this is the first iteration where a max MaxEntScan hasn't been set yet or the current value is larger than the current maximum
		if [ "$max_maxentscan" == "" ] || (( $(echo $current_maxentscan'>'$max_maxentscan | bc -l) )); then
			max_maxentscan=$current_maxentscan
			max_maxentscan_sequence=$current_splice_site
			max_maxentscan_coordinates=$current_coordinate_range
		fi
	done
	
	# For the maximal MaxEntScan window sequence, calculate MaxEntScan for the germline sequence
	current_splice_site_germline=$(samtools faidx $reference_genome $max_maxentscan_coordinates | grep -v '^>')
	current_maxentscan_germline=$(echo $current_splice_site_germline | perl MaxEntScan/score5.pl - | cut -f 2)
	
	echo "$line"$'\t'"$max_maxentscan"$'\t'"$max_maxentscan_sequence"$'\t'"$max_maxentscan_coordinates"$'\t'"$current_maxentscan_germline"$'\t'"$current_splice_site_germline" >> $out_dir/$prefix.introme.annotated.tsv
	
done

echo $(date +%x_%r) 'MaxEntScan calculation complete'

##################################

exit
