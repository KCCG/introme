#!/usr/bin/perl

use strict;
use warnings;
use Tie::File; # File -> array for parsing without loading whole thing into RAM

####################################

my $input_file; # Stores the input filename from the argument passed to the script

if (scalar(@ARGV) != 1) {
	print "FATAL ERROR: arguments must be supplied as 1) input file path.\n";
	exit;
} else {
	$input_file = $ARGV[0];
}

####################################

# Check that the file exists
-e $input_file or die "File \"$input_file\" does not exist.\n";

my @input_file_lines;

# Load the input file into an array (each element is a line but the whole thing is not loaded into memory)
tie @input_file_lines, 'Tie::File', $input_file or die "Cannot index input file.\n";

####################################

my $outputfile = $input_file.".vcf";

open(OUTPUT, ">$outputfile") or die "Cannot open file \"$outputfile\" to write to.\n\n";

####################################

for (my $i = 0; $i < scalar(@input_file_lines); $i++) {

	# Ignore header lines
	if ($input_file_lines[$i] =~ /^chromosome.*/) {
		print OUTPUT "##fileformat=VCFv4.1\n";
		print OUTPUT "##INFO=<ID=strand,Number=1,Type=String,Description=\"strand\">\n";
		print OUTPUT "##INFO=<ID=exon_id,Number=1,Type=String,Description=\"exon_id\">\n";
		print OUTPUT "##INFO=<ID=exon_number,Number=1,Type=Integer,Description=\"exon_number\">\n";
		print OUTPUT "##INFO=<ID=exon_5prime,Number=1,Type=String,Description=\"exon_5prime\">\n";
		print OUTPUT "##INFO=<ID=to_5prime,Number=1,Type=Integer,Description=\"to_5prime\">\n";
		print OUTPUT "##INFO=<ID=to_3prime,Number=1,Type=Integer,Description=\"to_3prime\">\n";
		print OUTPUT "##INFO=<ID=seq_motif,Number=1,Type=String,Description=\"seq_motif\">\n";
		print OUTPUT "##INFO=<ID=branchpoint_nt,Number=1,Type=String,Description=\"branchpoint_nt\">\n";
		print OUTPUT "##INFO=<ID=branchpoint_prob,Number=1,Type=Float,Description=\"branchpoint_prob\">\n";
		print OUTPUT "##INFO=<ID=U2_binding_energy,Number=1,Type=Float,Description=\"U2_binding_energy\">\n";
		print OUTPUT "##INFO=<ID=in_testtrain,Number=1,Type=Integer,Description=\"in_testtrain\">\n";
		print OUTPUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
		
		next;
	}
	
	# Split the row by comma characters
	my @split_line = split(/,/, $input_file_lines[$i]);
	
	if (scalar(@split_line) == 0) {
		print "WARNING: Found a line in the input file that does not contain comma-separated values. Line number: ".($i + 1)." Line contents: ".$input_file_lines[$i]."\n";
		
		next;
	} elsif (scalar(@split_line) != 14) {
		print "WARNING: Found a line that doesn't contain 14 columns of information. Line number: ".($i + 1)." Line contents: ".$input_file_lines[$i]." Number of columns detected: ".scalar(@split_line)."\n";
		
		next;
	}
	
	my $chromosome = $split_line[0];
	my $start = $split_line[1];
	my $end = $split_line[2];
	my $strand = $split_line[3];
	my $exon_id = $split_line[4];
	my $exon_number = $split_line[5];
	my $exon_5prime = $split_line[6];
	my $to_5prime = $split_line[7];
	my $to_3prime = $split_line[8];
	my $seq_motif = $split_line[9];
	my $branchpoint_nt = $split_line[10];
	my $branchpoint_prob = $split_line[11];
	my $U2_binding_energy = $split_line[12];
	my $in_testtrain = $split_line[13];
	
	my $ref = "N";
	my @alts;
	
	# Remove chr from the start of chromosome names
	#$chromosome =~ s/^chr//;
	
	# If the branchpoint is on the negative strand, the reference allele needs to be reversed to match the reference genome
	if ($strand eq "-") {
		if ($branchpoint_nt eq "A") {
			$ref = "T";
		} elsif ($branchpoint_nt eq "T") {
			$ref = "A";
		} elsif ($branchpoint_nt eq "C") {
			$ref = "G";
		} elsif ($branchpoint_nt eq "G") {
			$ref = "C";
		}
	} else {
		$ref = $branchpoint_nt;
	}
	
	# Determine all possible alt alleles for the current branchpoint allele
	if ($ref eq "A") {
		push(@alts, "T", "C", "G");
	} elsif ($ref eq "T") {
		push(@alts, "A", "C", "G");
	} elsif ($ref eq "C") {
		push(@alts, "A", "G", "T");
	} elsif ($ref eq "G") {
		push(@alts, "A", "T", "C");
	} else {
		die "FATAL ERRROR: could not determine reference allele or reference allele is not A, C, T, G. Line: ".$input_file_lines[$i];
	}
	
	foreach my $alt (@alts) {
		print OUTPUT $chromosome."\t"; # CHROM
		print OUTPUT $start."\t"; # POS
		print OUTPUT ".\t"; # ID	
		print OUTPUT $ref."\t"; # REF
		print OUTPUT $alt."\t"; # ALT
		print OUTPUT "1\t"; # QUAL
		print OUTPUT "PASS\t"; # FILTER
		print OUTPUT "strand=".$strand.";exon_id=".$exon_id.";exon_number=".$exon_number.";exon_5prime=".$exon_5prime.";to_5prime=".$to_5prime.";to_3prime=".$to_3prime.";seq_motif=".$seq_motif.";branchpoint_nt=".$branchpoint_nt.";branchpoint_prob=".$branchpoint_prob.";U2_binding_energy=".$U2_binding_energy.";in_testtrain=".$in_testtrain."\n"; # INFO
	}
		
	if ($i =~ /0000$/) {
		print "[".localtime()."] Processed: ".$i." lines from a total of ".(scalar(@input_file_lines)-1).".\n";
	}
}

close OUTPUT;

print "Done\n\n";

exit;