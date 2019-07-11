#!/usr/bin/perl

# Requirements:	Languages: 	perl
# Usage:			perl create_pool.pl [LIST_OF_GENOMES-FILE] [NUMBER_OF_GENOMES_FROM_FILE_THAT_SHOULD_BE_POOLED_TOGETHER] [NUMBER_OF_BASES_THAT_SHOULD_BE_POOLED_TOGETHER_IN_TOTAL] [PREFIX]
# Example:		perl create_pool.pl example1_list_of_plasmids.txt 3 10000000 example1_plasmid
use strict;
use warnings;
use List::Util qw/shuffle/;
use Cwd 'abs_path';

my ($i,$number_of_genomes,$total_size,$reads_per_genome,$genome,$counter,$line1,$read);		# initiate variables
my (@all_genomes,@chosen_genomes,@all_reads,@all_reads_shuffled);

$number_of_genomes=$ARGV[1];	# save genomes to shuffle	
$total_size=$ARGV[2];		# save total bases
$reads_per_genome=int(($total_size/$number_of_genomes)/8370);	#calculate pulled reads per genome

open(STATS,">","$ARGV[3]_stats.txt");						# open output file for stats
print STATS "List of All Genomes:\t\t$ARGV[0]\nIDs And Pathways:\t\t$ARGV[3]_ids_and_pathways.txt\nRead Pool:\t\t\t$ARGV[3]_read_pool.txt\nNumber Of Chosen Genomes:\t$ARGV[1]\nTotal Size:\t\t\t$ARGV[2]\nNumber Of Reads Per Genome:\t$reads_per_genome\n";	# write info in stats file
print "\n\t- - -   Pulling $ARGV[1] random genomes outof list ($ARGV[0])   - - -\n\n";

open(INONE,"$ARGV[0]") or die "File $ARGV[0] does not exist";	# open list of genome files
@all_genomes=shuffle<INONE>;								# shuffle list into an array
close(INONE);

open(OUTTWO,">","$ARGV[3]_ids_and_pathways.txt");				# open first output file
for($i=0;$i<$number_of_genomes;$i++){						# for the first X genomes
	chomp $all_genomes[$i];
	$all_genomes[$i]=~ m/(.+)(\.fna)/g;					# filter ID outof file name 
	push(@chosen_genomes,$1);							# push ID into new array
	my $root1 = abs_path("Sim_Pipeline/$1_l8500_c150/$1_l8500_c150-filtered_R1.fastq");	# save absolute pathway
	my $root2 = abs_path("Sim_Pipeline/$1_l8500_c150/$1_l8500_c150-filtered_R2.fastq");	# save absolute pathway
	print OUTTWO "$1\t$root1\t$root2\n"; 										# print ID and illumina pathways into first output file
}										
close OUTTWO;
print "\n\t- - -   Pulling $reads_per_genome random reads per genome outof PacBio read files   - - -\n";

open(OUTTHREE,">","read_pool.tmp");				# open temporary output file
foreach $genome (@chosen_genomes){					# parse IDs of all chosen genomes
	print "\t\tProcessing $genome...\n";
	$counter=0;								# set counter (back) to 0
	undef(@all_reads);							# undefine array from last iteration
	undef(@all_reads_shuffled);					# undefine array from last iteration
print "/gpfs/project/dilthey/projects/multiplexer_pool/Ultraplexing_example/Sim_Pipeline/".$genome."_l8500_c150/".$genome."_l8500_c150-filtered_complete.fastq\n";
	open(INTWO,"Sim_Pipeline/".$genome."_l8500_c150/".$genome."_l8500_c150-filtered_complete.fastq") or die "Read-File does not exist";	# open PacBio read file
	$read="@".$genome."_";						# save ID	in variable									
	while($line1=<INTWO>){						# parse read file
		chomp($line1);
		$counter++;
		if(($counter%4)==0){					# every 4th line (4 lines for one read each)
			$read.="$line1";					# add line to variable
			push(@all_reads,$read);				# push variable into array
			$read="@".$genome."_";				# reset variable to ID
		}
		else{								# line 1-3
			$read.="$line1>>>>>NEWLINE>>>>>";		# add line to variable	# add 10x"-" (for later replacement)
		}
	}
	@all_reads_shuffled=shuffle(@all_reads);			# shuffle array (all reads of one genome) into new array
	if ($reads_per_genome<=scalar(@all_reads_shuffled)){
		for($i=0;$i<$reads_per_genome;$i++){			# for a certain number of reads
			print OUTTHREE "$all_reads_shuffled[$i]\n";	# print read into temporary output file
		}
	}
	else{ die "\n\n\t\tERROR:\tMore reads claimed ($reads_per_genome) than available (".scalar(@all_reads_shuffled).").\n\t\t\tProcess will be aborted.\n\t\t\tPlease lower the total number of basepairs that should be taken into the pool )\n";
	}
}
close(INTWO);
print "\n\t- - -   Final shuffling and cleanup   - - -\n\n";
system("shuf read_pool.tmp > read_pool_shuffled.tmp");										# shuffle temporary output file
system("perl -p -e 's/>>>>>NEWLINE>>>>>/\n/g;' read_pool_shuffled.tmp > $ARGV[3]_read_pool.fastq");	# replace 10x"-" with a newline "\n"

system("rm read_pool.tmp");														# remove temporary files
system("rm read_pool_shuffled.tmp");												# remove temporary files
undef(@all_genomes);															# empty array
undef(@chosen_genomes);															# empty array
undef(@all_reads);																# empty array
undef(@all_reads_shuffled);														# empty array
										


