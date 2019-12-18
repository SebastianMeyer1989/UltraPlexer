#!/usr/bin/perl
# FOR REAL DATA (data without original genome)
# Usage:	perl create_kmer_based_fastq_for_real_data.pl [CALLED_KMERS] [USED_READPOOL] [PREFIX]

#use strict;
use warnings;
use List::MoreUtils qw/uniq/;

if(scalar(@ARGV)<3){die "Please specify parameters on command line (see script)."}	# abort script if to few parameters

my ($line1,$line2,$prefix,$file);		# initiate variables
my (%called_kmers,@files);

$prefix=$ARGV[2];										# prefix for the run
system("rm -rf $prefix*predicted_reads.fastq");				# remove old files with identical prefix
system("rm -rf $prefix*true_reads.fastq");					# remove old files with identical prefix

open(INONE,"$ARGV[0]") or die "File $ARGV[0] does not exist";	# Open File: [prefix].called_kmers / [prefix].called_random
while($line1=<INONE>){									# parse calling table
	chomp($line1);
	if($line1=~m/^(.+)\t(.+)\t(.+)$/g){					# 1st Column = Original, 2nd Column = Prediction
		$called_kmers{"\@$1"}="$2";						# save predicted id under original id
	}
}
close(INONE);

open(INTWO,"$ARGV[1]") or die "File $ARGV[1] does not exist";	# Open File: [prefix]_readpool.fastq
while($line2=<INTWO>){									# parse read pool
	chomp($line2);
	if($line2=~m/^(@.+?) .+/){	# if header has empty spaces
		$head=$1;				# ... save first part of header
	}
	elsif($line2=~m/^@/){		# else
		$head=$line2;			# save header
	}
	if(($line2=~m/^@/) && (exists($called_kmers{$head}))){					# if line is a header and read exists in called kmers
		if($called_kmers{$head} ne 'NA')
		{
			push(@files, "$called_kmers{$head}");							# save id
			open(TMP1,">>","$prefix-$called_kmers{$head}-predicted_reads.fastq");	# Output File: Different fastq files, Predictions, 2nd Column

			print TMP1 "$line2\n";		# write next 1/4 lines in prediction file (equal to 1 read)

			$line2=<INTWO>;			# read next line
			print TMP1 "$line2";		# write next 2/4 lines in prediction file (equal to 1 read)

			$line2=<INTWO>;			# read next line
			print TMP1 "$line2";		# write next 3/4 lines in prediction file (equal to 1 read)

			$line2=<INTWO>;			# read next line
			print TMP1 "$line2";		# write next 4/4 lines in prediction file (equal to 1 read)

			close(TMP1);	
		}

	}
	elsif(($line2=~m/^@/) && (!exists($called_kmers{$head}))){		# if line is a header and read does NOT exist in called kmers
		print "Read -- $head -- not in list of calles kmers.\n";	
	}
	else{}
}

