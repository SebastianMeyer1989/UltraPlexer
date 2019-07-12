#!/usr/bin/env perl

# Requirements:	Languages: 	perl
# Usage:			perl parse_calling_tbl.pl [CALLING_TABLE]
# Example:		perl simulation_pipeline.pl example1_plasmid.classification_k19.called_kmers
#

use strict;
use warnings;

my ($line1,$count1,$count2,%ratiocount,$key);
$count2=0;
open(OUTONE,">","$ARGV[0].tmp");
open(IN,"$ARGV[0]");					# open calling table
while($line1=<IN>){
	chomp($line1);
	if($line1=~m/^(.+)_@.+\t(.+)\t.+/){
		if($1 eq $2){					# if assigned correctly
			$ratiocount{"$1\ttrue"}++;	# save in hash
			$count1++;
print OUTONE "$1\t$2\tcorrect\n";	
		}
		else{						# if not assigned correctly
			$ratiocount{"$1\tfalse"}++;	# save in hash
			$count2++;
print OUTONE "$1\t$2\tfalse\n";
		}
	}
}
close(IN);

print OUTONE "$ARGV[0]\tCorrect_Reads:\t$count1\tFalse_Reads:\t$count2\tRatio_Correct_Reads:\t".$count1/($count1+$count2)."\n";
system("less $ARGV[0].tmp | sort -k3 | uniq -c | sort -k3  > $ARGV[0].stats");
close(OUTONE);

system("rm $ARGV[0].tmp");

open(OUTTHREE,">","$ARGV[0].stats2");		# create stats file
print OUTTHREE "$ARGV[0] Summary\t\tCorrect_Reads:\t$count1\tFalse_Reads:\t$count2\tRatio_Correct_Reads:\t".$count1/($count1+$count2)."\n";		# write summary to stats file
foreach $key(sort(keys(%ratiocount))){
	if($key=~m/.+false/){							
		print OUTTHREE "$key: $ratiocount{$key}\t\t";	# first print false qssignment,
	}
	else{
		print OUTTHREE "$key: $ratiocount{$key}\n";		# ... then correct assignment
	}
}
