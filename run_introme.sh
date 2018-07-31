#!/bin/bash

# Introme Bash Script
#
# Developers:
# Velimir Gayevskiy (vel@vel.nz)
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)
#





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
CADD_Phred='>=0' # Minimum CADD score (phred-scaled)
min_QUAL='>=200' # The QUAL VCF field
min_DP='>=20' # The sample with the highest depth must have a equal or larger depth than this value
min_AD='>=5' # The sample with the highest number of alternate reads must have a equal or larger number than this value

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
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

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
while [[ $i -le $(( ${#vcfsamples[@]} - 1 )) ]]; do 
	# Code to not apply genotype filtering to unaffected samples, will implement properly soon...
	#stringinarray "${vcfsamples[$i]}" "${affected_samples[@]}"
	
	#if [ $? != 0 ]; then
	#	(( i++ ))
		
	#	continue
	#fi
	
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
	(( i++ ))
done

# Cut the 4 first characters off the filter string (the " && ")
genotype_filters=$(echo $genotype_filters | cut -c 4-)

bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"$genotype_filters" $out_dir/$prefix.subset.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.vcf.gz
tabix $out_dir/$prefix.subset.inheritancefilter.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.inheritancefilter.vcf.gz | grep -v '^#' | wc -l) 'variants after familial filter'
echo $(date +%x_%r) 'Familial filtering complete'

##################################
# STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream

echo $(date +%x_%r) 'Beginning annotation'

vcfanno -p $(getconf _NPROCESSORS_ONLN) conf.toml $out_dir/$prefix.subset.inheritancefilter.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
tabix -p vcf $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'

##################################
# STEP 4: Hard filtering the subsetted, annotated VCF

echo $(date +%x_%r) 'Beginning filtering'

bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && TYPE='snp' && (QUAL$min_QUAL || QUAL='.') && MAX(FORMAT/DP[*])$min_DP && MAX(FORMAT/AD[*:1])$min_AD && (MGRB_AF$MGRB_AF || MGRB_AF='.') && (gnomAD_PM_AF$gnomad_popmax_AF || gnomAD_PM_AF='.') && (Branchpointer_Branchpoint_Prob!='.' || (Branchpointer_Branchpoint_Prob='.' && (CADD_Phred$CADD_Phred || CADD_Phred='.')))" $out_dir/$prefix.subset.inheritancefilter.annotated.vcf.gz | bgzip > $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz
tabix -p vcf $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz

echo $(gzip -d -c $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Filtering complete'

##################################
# STEP 5: Output final list of variants as a spreadsheet-friendly TSV

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/Gene_Symbol\t%INFO/Gene_Strand\t%INFO/CADD_Phred\t%INFO/MGRB_AF\t%INFO/gnomAD_PM_AF\t%INFO/CSQ\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/Branchpointer_Branchpoint_Prob\t%INFO/Branchpointer_U2_Binding_Energy[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/$prefix.subset.inheritancefilter.annotated.filtered.vcf.gz > $out_dir/$prefix.introme.tsv

# Remove the '[<number]' column numbers added by bcftools query to column names
grep '^#' $out_dir/$prefix.introme.tsv | sed "s/\[[0-9]*\]//g" > $out_dir/$prefix.introme.tsv.header
grep -v '^#' $out_dir/$prefix.introme.tsv >> $out_dir/$prefix.introme.tsv.header
mv $out_dir/$prefix.introme.tsv.header $out_dir/$prefix.introme.tsv

# Sort by chromosome and coordinate
sort -k1,1n -k2,2n $out_dir/$prefix.introme.tsv > $out_dir/$prefix.introme.sorted.tsv
mv $out_dir/$prefix.introme.sorted.tsv $out_dir/$prefix.introme.tsv

##################################
# STEP 6: Calculate the maximum 3' and 5' MaxEntScan value for sliding windows around each variant

echo $(date +%x_%r) 'Beginning MaxEntScan calculation, number of variants to process: '$(cat $out_dir/$prefix.introme.tsv | wc -l)

# Number of lines processed for displaying progress to the user
num_lines_processed=0

# Go through each line of the results TSV
cat $out_dir/$prefix.introme.tsv | \

while read line; do
	# If the current number of running processes is below a certain number, wait to process further lines
	# Processing each line uses 107 processes on one core (currently, this can change!) so use this as an indicator of how many cores this will simultaneously use
	while [[ $(jobs -r | wc -l) -ge 450 ]]; do
		sleep 1
	done
	
	# Iterate the number of lines processed
	((num_lines_processed++))
	
	# Report progress for every 10 variants processed
	if [[ $(echo $num_lines_processed | grep '0$') ]]; then
		echo $(date +%x_%r) 'Processed '$num_lines_processed' variants'
	fi
	
	# Create a block for processing and calculating MaxEntScan scores which is pushed into the background at the end
	{
		# If the line is the header line
		if [[ ${line:0:1} == "#" ]]; then
			# Add to the header line
			echo "$line"$'\t'"MaxEntScan_Consequence"$'\t'"5'_max_MaxEntScan_value"$'\t'"5'_max_MaxEntScan_reference_value"$'\t'"3'_max_MaxEntScan_value"$'\t'"3'_max_MaxEntScan_reference_value"$'\t'"5'_max_MaxEntScan_direction"$'\t'"3'_max_MaxEntScan_direction"$'\t'"5'_max_MaxEntScan_coordinates"$'\t'"3'_max_MaxEntScan_coordinates"$'\t'"5'_max_MaxEntScan_sequence"$'\t'"5'_max_MaxEntScan_reference_sequence"$'\t'"3'_max_MaxEntScan_sequence"$'\t'"3'_max_MaxEntScan_reference_sequence" > $out_dir/$prefix.introme.annotated.tsv
		
			# Don't do any further processing on this line
			continue
		fi
	
		# Current variant information
		current_chr=$(echo "$line" | cut -f 1)
		current_pos=$(echo "$line" | cut -f 2)
		current_alt=$(echo "$line" | cut -f 4)
		current_strand=$(echo "$line" | cut -f 7)
	
		# Calculate the maximum range we are interested in (corresponding to the 23-base windows)
		current_max_coordinate_range="$current_chr:"$(( $current_pos - 22 ))"-"$(( $current_pos + 22 ))
	
		# Fetch the reference genome sequence once and subset it for windows
		current_max_reference_context=$(samtools faidx $reference_genome $current_max_coordinate_range | grep -v '^>')
		
		# Determine the direction in which to look for putative splice sites based on the direction of the genes with which the variant overlaps
		# If the variant overlaps with both a + and - strand gene
		if [[ $(echo $current_strand | grep '-') ]] && [[ $(echo $current_strand | grep '+') ]]; then
			current_strand_direction_calculation='unknown'
		# If the variant only overlaps with - strand gene(s)
		elif [[ $(echo $current_strand | grep '-') ]]; then
			current_strand_direction_calculation='reverse'
		# If the variant only overlaps with + strand gene(s)
		elif [[ $(echo $current_strand | grep '+') ]]; then
			current_strand_direction_calculation='forward'
		# If the variant does not have a strand direction associated (it probably doesn't overlap with a gene)
		else
			current_strand_direction_calculation='unknown'
		fi
		
		##################################
	
		# 5' MaxEntScan score calculation
	
		max_maxentscan_5=''
		max_maxentscan_sequence_5=''
		max_maxentscan_coordinates_5=''
		max_maxentscan_difference_5=''
		max_maxentscan_consequence_5=''
		max_maxentscan_type_5=''
	
		# $i is the offset for determining the coordinate range
		for i in {0..8}; do
			# Determine the putative splice site coordinate range
			current_coordinate_range="$current_chr:"$(( $current_pos - 8 + $i ))"-"$(( $current_pos + $i ))
			
			# Extract the current putative splice site sequence as in the reference genome and mutate the correct base corresponding to the variant
			current_splice_site_forward=${current_max_reference_context:(( 14 + $i )):(( 8 - $i ))}$current_alt${current_max_reference_context:23:$i} # The numbers here are offsetting to ignore the full 23 base window and focus on the 9 base window
			
			# If the variant gene is not on the reverse strand, calculate the MaxEntScan value for the splice site being on the forward strand
			if [[ $current_strand_direction_calculation != "reverse" ]]; then
				current_maxentscan_forward=$(echo $current_splice_site_forward | perl MaxEntScan/score5.pl - | cut -f 2)
			fi
			
			# If the variant gene is not on the forward strand, calculate reverse values too
			if [[ $current_strand_direction_calculation != "forward" ]]; then
				# Determine the reverse complement putative splice site sequence
				current_splice_site_reverse=$(echo $current_splice_site_forward | tr "[ATCG]" "[TAGC]" | rev)
			
				# Calculate the MaxEntScan value for the splice site being on the reverse strand
				current_maxentscan_reverse=$(echo $current_splice_site_reverse | perl MaxEntScan/score5.pl - | cut -f 2)
			fi
			
			# If the variant gene is on the forward strand, or the strand is unknown but the forward value is higher than the reverse value
			if [[ $current_strand_direction_calculation == "forward" ]] || ([[ $current_strand_direction_calculation == "unknown" ]] && (( $(echo $current_maxentscan_forward' > '$current_maxentscan_reverse | bc -l) ))); then
				max_maxentscan_5=$current_maxentscan_forward
				max_maxentscan_sequence_5=$current_splice_site_forward
				max_maxentscan_coordinates_5=$current_coordinate_range
				max_maxentscan_type_5="forward"
			else
				max_maxentscan_5=$current_maxentscan_reverse
				max_maxentscan_sequence_5=$current_splice_site_reverse
				max_maxentscan_coordinates_5=$current_coordinate_range
				max_maxentscan_type_5="reverse"
			fi
		done
		
		# Fetch the reference genome sequence for the window coordinates
		max_maxentscan_sequence_5_reference=$(samtools faidx $reference_genome $max_maxentscan_coordinates_5 | grep -v '^>')
	
		# Adjust the reference genome sequence if the most damaging putative splice site is on the reverse strand
		if [[ $max_maxentscan_type_5 == "reverse" ]]; then
			max_maxentscan_sequence_5_reference=$(echo $max_maxentscan_sequence_5_reference | tr "[ATCG]" "[TAGC]" | rev)
		fi
		
		# Calculate the reference genome score
		max_maxentscan_5_reference=$(echo $max_maxentscan_sequence_5_reference | perl MaxEntScan/score5.pl - | cut -f 2)
	
		# Calculate the absolute difference between the reference and variant scores (i.e. take negative numbers into account)
		max_maxentscan_difference_5=$(echo $max_maxentscan_5" - "$max_maxentscan_5_reference | bc -l | sed 's/-//') # Force the difference to be positive so it's an absolute difference
	
		# Determine the potential MaxEntScan consequence (i.e. likelihood of a new splice site to replace the canonical one)
		# If the score is below the reference score or below zero
		if (( $(echo $max_maxentscan_5' <= '$max_maxentscan_5_reference | bc -l) )) || (( $(echo $max_maxentscan_5' < 0' | bc -l) )); then
			max_maxentscan_consequence_5='NONE'
		# If the score is 0-4 and the difference to the reference score is >=4
		elif (( $(echo $max_maxentscan_5' < 4' | bc -l) )) && (( $(echo $max_maxentscan_difference_5' >= 4' | bc -l) )); then
			max_maxentscan_consequence_5='LOW'
		# If the score is 4-10 and the difference to the reference score is >=4
		elif (( $(echo $max_maxentscan_5' < 10' | bc -l) )) && (( $(echo $max_maxentscan_difference_5' >= 4' | bc -l) )); then
			max_maxentscan_consequence_5='MED'
		# If the score is greater than 10 and the difference to the reference score is >= 6
		elif (( $(echo $max_maxentscan_5' >= 10' | bc -l) )) && (( $(echo $max_maxentscan_difference_5' >= 6' | bc -l) )); then
			max_maxentscan_consequence_5='HIGH'
		fi
	
		##################################
	
		# 3' MaxEntScan score calculation
	
		max_maxentscan_3=''
		max_maxentscan_sequence_3=''
		max_maxentscan_coordinates_3=''
		max_maxentscan_difference_3=''
		max_maxentscan_consequence_3=''
		max_maxentscan_type_3=''
	
		# $i is the offset for determining the coordinate range
		for i in {0..22}; do
			# Determine the putative splice site coordinate range
			current_coordinate_range="$current_chr:"$(( $current_pos - 22 + $i ))"-"$(( $current_pos + $i ))
		
			# Extract the current putative splice site sequence as in the reference genome and mutate the correct base corresponding to the variant
			current_splice_site_forward=${current_max_reference_context:(( 0 + $i )):(( 22 - $i ))}$current_alt${current_max_reference_context:23:$i}
			
			# If the variant gene is not on the reverse strand, calculate the MaxEntScan value for the splice site being on the forward strand
			if [[ $current_strand_direction_calculation != "reverse" ]]; then
				current_maxentscan_forward=$(echo $current_splice_site_forward | perl MaxEntScan/score3.pl - | cut -f 2)
			fi
			
			# If the variant gene is not on the forward strand, calculate reverse values too
			if [[ $current_strand_direction_calculation != "forward" ]]; then
				# Determine the reverse complement putative splice site sequence
				current_splice_site_reverse=$(echo $current_splice_site_forward | tr "[ATCG]" "[TAGC]" | rev)
			
				# Calculate the MaxEntScan value for the splice site being on the reverse strand
				current_maxentscan_reverse=$(echo $current_splice_site_reverse | perl MaxEntScan/score3.pl - | cut -f 2)
			fi
			
			# If the variant gene is on the forward strand, or the strand is unknown but the forward value is higher than the reverse value
			if [[ $current_strand_direction_calculation == "forward" ]] || ([[ $current_strand_direction_calculation == "unknown" ]] && (( $(echo $current_maxentscan_forward' > '$current_maxentscan_reverse | bc -l) ))); then
				max_maxentscan_3=$current_maxentscan_forward
				max_maxentscan_sequence_3=$current_splice_site_forward
				max_maxentscan_coordinates_3=$current_coordinate_range
				max_maxentscan_type_3="forward"
			else
				max_maxentscan_3=$current_maxentscan_reverse
				max_maxentscan_sequence_3=$current_splice_site_reverse
				max_maxentscan_coordinates_3=$current_coordinate_range
				max_maxentscan_type_3="reverse"
			fi
		done
	
		# For the maximal MaxEntScan window sequence, calculate MaxEntScan for the reference sequence
		max_maxentscan_sequence_3_reference=$(samtools faidx $reference_genome $max_maxentscan_coordinates_3 | grep -v '^>')
		
		# Adjust the reference genome sequence if the most damaging putative splice site is on the reverse strand
		if [[ $max_maxentscan_type_3 == "reverse" ]]; then
			max_maxentscan_sequence_3_reference=$(echo $max_maxentscan_sequence_3_reference | tr "[ATCG]" "[TAGC]" | rev)
		fi
		
		# Calculate the reference genome score
		max_maxentscan_3_reference=$(echo $max_maxentscan_sequence_3_reference | perl MaxEntScan/score3.pl - | cut -f 2)
	
		# Calculate the absolute difference between the reference and variant scores (i.e. take negative numbers into account)
		max_maxentscan_difference_3=$(echo $max_maxentscan_3" - "$max_maxentscan_3_reference | bc -l | sed 's/-//') # Force the difference to be positive so it's an absolute difference
	
		# Determine the potential MaxEntScan consequence (i.e. likelihood of a new splice site to replace the canonical one)
		# If the score is below the reference score or below zero
		if (( $(echo $max_maxentscan_3' <= '$max_maxentscan_3_reference | bc -l) )) || (( $(echo $max_maxentscan_3' < 0' | bc -l) )); then
			max_maxentscan_consequence_3='NONE'
		# If the score is 0-4 and the difference to the reference score is >=4
		elif (( $(echo $max_maxentscan_3' < 4' | bc -l) )) && (( $(echo $max_maxentscan_difference_3' >= 4' | bc -l) )); then
			max_maxentscan_consequence_3='LOW'
		# If the score is 4-10 and the difference to the reference score is >=4
		elif (( $(echo $max_maxentscan_3' < 10' | bc -l) )) && (( $(echo $max_maxentscan_difference_3' >= 4' | bc -l) )); then
			max_maxentscan_consequence_3='MED'
		# If the score is greater than 10 and the difference to the reference score is >= 6
		elif (( $(echo $max_maxentscan_3' >= 10' | bc -l) )) && (( $(echo $max_maxentscan_difference_3' >= 6' | bc -l) )); then
			max_maxentscan_consequence_3='HIGH'
		fi
	
		##################################
	
		# Determine the maximum MaxEntScan consequence
	
		max_maxentscan_consequence=''
	
		if [[ $max_maxentscan_consequence_5 == 'HIGH' ]] || [[ $max_maxentscan_consequence_3 == 'HIGH' ]]; then
			max_maxentscan_consequence='HIGH'
		elif [[ $max_maxentscan_consequence_5 == 'MED' ]] || [[ $max_maxentscan_consequence_3 == 'MED' ]]; then
			max_maxentscan_consequence='MED'
		elif [[ $max_maxentscan_consequence_5 == 'LOW' ]] || [[ $max_maxentscan_consequence_3 == 'LOW' ]]; then
			max_maxentscan_consequence='LOW'
		else
			max_maxentscan_consequence='NONE'
		fi
	
		##################################
		
		# Scatter the output variants to separate temporary files
		
		# Create an empty temporary file
		tmpfile=$(mktemp /tmp/$prefix.introme.annotated.tsv.XXXXX)
		
		# Output the result to the temporary file
		echo "$line"$'\t'"$max_maxentscan_consequence"$'\t'"$max_maxentscan_5"$'\t'"$max_maxentscan_5_reference"$'\t'"$max_maxentscan_3"$'\t'"$max_maxentscan_3_reference"$'\t'"$max_maxentscan_type_5"$'\t'"$max_maxentscan_type_3"$'\t'"$max_maxentscan_coordinates_5"$'\t'"$max_maxentscan_coordinates_3"$'\t'"$max_maxentscan_sequence_5"$'\t'"$max_maxentscan_sequence_5_reference"$'\t'"$max_maxentscan_sequence_3"$'\t'"$max_maxentscan_sequence_3_reference" >> $tmpfile
	} &	
done

# Gather all temporary files into the results file
for file in /tmp/$prefix.introme.annotated.tsv.*; do 
	cat "$file"
done >> $out_dir/$prefix.introme.annotated.tsv

echo $(date +%x_%r) 'MaxEntScan calculation complete'

# Sort by chromosome and coordinate
sort -k1,1n -k2,2n $out_dir/$prefix.introme.annotated.tsv > $out_dir/$prefix.introme.annotated.sorted.tsv
mv $out_dir/$prefix.introme.annotated.sorted.tsv $out_dir/$prefix.introme.annotated.tsv

##################################

exit
