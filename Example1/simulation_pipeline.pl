#!/usr/bin/env perl

# Requirements:	Languages: 	perl
# Requirements:	Programmes:	pbsim, wgsim,
# Usage:			perl simulation_pipeline.pl [LIST_OF_GENOMES-FILE] [NANOPORE_FRAGMENT_LEGNTH] [NANOPORE_COVERAGE]
# Example:		perl simulation_pipeline.pl example1_list_of_plasmids.txt 8500 150
# Important:		- This script needs to be in the same folder as the genome-files
#				- Files in list must have the following format: [NAME].[ENDING]     (Examples: MRSA252.fna, LGA251.fasta, wrongfile.txt, s_aureus.fa)
#				- Pathways for simulations need to be changed in the following three lines:
my $pbsimpath="/gpfs/project/dilthey/software/PBSIM-PacBio-Simulator/src/pbsim";		# path of pbsim (long-read simulation
my $modelpath="/gpfs/project/dilthey/software/PBSIM-PacBio-Simulator/data/model_qc_clr";	# qc-model for pbsim
my $wgsimpath="/gpfs/project/dilthey/software/wgsim/wgsim";							# path of wgsim (short-read simulation)

use strict;
use warnings;
use Time::localtime;


################################################################################################################
##########  Initiate Variables And Co.  ########################################################################
################################################################################################################

my ($tm,$d,$mon,$y,$h,$min,$sec,$i);																					# initiate time variables
my ($line1,$line2,$line3,$line4,$line5,$line6,$line7,$string1,$string2,$r1line,$r1in,$r2in,$count2,$count,$count3,$count4,$count5,$max,%sizehash1,%sizehash2);		# initiate variables
my ($break,$pre,$r1,$r2,$lr,$maf,$sub,$spades,$racon,$pilon,$header,$conca,$mafline,$pbin,$ident,$cover,$ident2,$cover2,$lsd,$p1,$p2,$header2,$counter,$key);		# initiate variables

system ("mkdir -p Sim_Pipeline");							# create output folder
open(OUTONE,">>","Sim_Pipeline/missing_files.log") or die "Can't create log-file"; 				# create file for missing lineages
open(OUTTWO,">>","Sim_Pipeline/simulation_parameters.log") or die "Can't create parameter-file"; 	# create file with simulation parameters

$tm=localtime;																# create date and time
($d,$mon,$y,$h,$min,$sec)=($tm->mday,$tm->mon,$tm->year,$tm->hour,$tm->min,$tm->sec);	# create date and time
$mon+=1;
$y+=1900;

$break= "____________________________________________________________________________________________________";	# save break line (for optical enhancement)
print OUTONE "\n\n$break\nMissing (or wrong spelled) files: $d.$mon.$y, $h:$min:$sec\n$break\n\n";				# header for missing files
print OUTTWO "\n$break\nPipeline Parameters: $d.$mon.$y, $h:$min:$sec\n$break\n\n";	 						# header for pipeline parameters
print OUTTWO "Output Pathway:\tSim_Pipeline/[FILE]/\n";
print OUTTWO "Simulated Reads:\t[FILE]_0001.fastq or .maf or .ref (PacBio), [FILE]_R1.fastq or R2 (Illumina)\n";
print OUTTWO "Filtered Reads:\t[FILE]-filtered_0001.fastq or .maf or .ref (PacBio), [FILE]-filtered_R1.fastq & [FILE]-filtered_R2.fastq or R2 (Illumina)\n";

die "Specify list of genomes as first argument" unless(-e $ARGV[0]);
die "Specify mean long-read fragment length as second argument" unless($ARGV[1] =~ /^\d+$/);
die "Specify long-read coverage as third argument" unless($ARGV[2] =~ /^\d+$/);

open(INONE,"$ARGV[0]") or die "Can't open $ARGV[0]";							# input list	 (names of gemome-files)
while($line1=<INONE>){							# parse input
	chomp($line1);
	next unless($line1);
	if(-e $line1){								# does file exist?
		die unless ($line1=~m/(.+)\.[a-zA-Z0-9]+$/);	# filter for part before "."
		$pre=$1."_l".$ARGV[1]."_c".$ARGV[2];		# define as prefix for output files
		$r1="_R1.fastq";						# save illumina-file-ending 1
		$r2="_R2.fastq";						# save illumina-file-ending 2
		if(-e "Sim_Pipeline/$pre")
		{
			system("rm -rf Sim_Pipeline/$pre") and die "Could not remove directory Sim_Pipeline/$pre";
		}
		system ("mkdir Sim_Pipeline/$pre");		# create seperate output folder for each lineage


################################################################################################################
##########  Parse Genome Files And Extend Sequences  ###########################################################
################################################################################################################

		open(OUTTHREE,">","Sim_Pipeline/$pre/$pre-extended.fna") or die "Can't create Sim_Pipeline/$pre/$pre-extended.fna";			# create file for extended genome
		open(INTWO, "$line1") or die "Can't open $line1";
		$conca="";												# clear concatenated genome variable
		$counter=0;

		while($line2=<INTWO>){										# parse genome file
			chomp($line2);
			if($line2=~m/^>(.+?) /){									# if header
				if($conca eq ""){									# ...and conca. variable is empty 
					print OUTTHREE "$line2-EXTENDED\n";				# print header (first of file)
					$header2=$1;									### header-ID speichern
					$counter++;									### header-Nummer speichern
				}
				else{											# ...and conca. variable is filled
					$sizehash1{$header2}=length($conca);				### größe zu header-ID speichern
					$sizehash2{$counter}=length($conca);				### größe zu header-Nummer speichern
					$string1 = substr($conca, 0, int(length($conca)/2));	# split conca. variable in half
					print OUTTHREE "$conca$string1\n";					# extend original genome and print it in file
					$conca="";									# clear concatenated genome variable
					print OUTTHREE "$line2-EXTENDED\n";				# print header (second to last of file)
					$header2=$1;									### header-ID speichern
					$counter++;										### header-Nummer speichern
				}
			}
			else {												# if not header
				$conca.=$line2;									# concatenate the genome sequence
			}
		}
		$sizehash1{$header2}=length($conca);							### größe zu header-ID speichern
		$sizehash2{$counter}=length($conca);							### größe zu header-Nummer speichern
		$string1 = substr($conca, 0, int(length($conca)/2));				# split conca. variable in half
		print OUTTHREE "$conca$string1\n";								# extend original genome and print it in file (last of file)
		$conca="";
		close(OUTTHREE);
		

################################################################################################################
##########  Do Simulations And Write Parameter Files  ##########################################################
################################################################################################################

#Simulation Options can be modified in this block


		$lsd=$ARGV[1]*0.76666;			# calculate length of standard deviation for pbsim
		my $max=$ARGV[1]*7.5;			# calculate maximum length for pbsim

		print "\n\n\n\n\t| | Start PacBio Simulation: $line1 | |\n\n\n\n";
		print OUTTWO "| | File: $pre | |\n";
		print OUTTWO "$pbsimpath --prefix Sim_Pipeline/$pre/$pre --depth $ARGV[2] --length-mean $ARGV[1] --length-sd $lsd --length-max $max --length-min 230  --accuracy-mean 0.88 --model_qc $modelpath Sim_Pipeline/$pre/$pre-extended.fna\n"; 	# write PacBio parameters to file
		system("$pbsimpath --prefix Sim_Pipeline/$pre/$pre --depth $ARGV[2] --length-mean $ARGV[1] --length-sd $lsd --length-max $max --length-min 230  --accuracy-mean 0.88 --model_qc $modelpath Sim_Pipeline/$pre/$pre-extended.fna");			# execute PacBio simulator		
		print "\n\n\n\n\t| | Start Illumina Simulation: $line1 | |\n\n\n\n";
		print OUTTWO "$wgsimpath -e 0.005 -1 150 -2 150 -d 278 -s 128 -r 0 -R 0 Sim_Pipeline/$pre/$pre-extended.fna Sim_Pipeline/$pre/$pre$r1 Sim_Pipeline/$pre/$pre$r2\n"; 	# write Illumina parameters to file
		system("$wgsimpath -e 0.005 -1 150 -2 150 -d 278 -s 128 -r 0 -R 0 Sim_Pipeline/$pre/$pre-extended.fna Sim_Pipeline/$pre/$pre$r1 Sim_Pipeline/$pre/$pre$r2");			# execute Illumina simulator


################################################################################################################
##########  Filter Read-Files And Erase Reads On Overhang  #####################################################
################################################################################################################

		foreach $key(keys(%sizehash2)){									### für alle sequencen, d.h. für alle long-read dateien (genom und plasmid)

			if(length($key)==1){
				$maf="_000".$key.".maf";									### save pb-file-ending maf
				$lr="_000".$key.".fastq";								### save pb-file-ending fastq
			}
			elsif(length($key)==2){
				$maf="_00".$key.".maf";									### save pb-file-ending maf
				$lr="_00".$key.".fastq";									### save pb-file-ending fastq
			}
			open(OUTFOUR,">>","Sim_Pipeline/$pre/$pre-filtered_complete.fastq") or die "Can't create Sim_Pipeline/$pre/$pre-filtered_complete.fastq";	### create file for filtered pb reads
			open(INFOUR,"Sim_Pipeline/$pre/$pre$maf") or die "Can't open Sim_Pipeline/$pre/$pre$maf";			# input maf file with length
			open(INFIVE,"Sim_Pipeline/$pre/$pre$lr") or die "Can't open Sim_Pipeline/$pre/$pre$lr";			# input pb reads for filtering
			<INFOUR>;
			while($mafline=<INFOUR>){							# parse maf file
				chomp($mafline);
				next unless($mafline);
				die unless($mafline=~m/\S+ \S+\s+(\d+)\s+(\d+)\s+/);
				if($mafline=~m/\S+ \S+\s+(\d+)\s+(\d+)\s+/){			# if reference line
					if($1>$sizehash2{$key} || $2>$sizehash2{$key}){	# and read starts after original / read is longer than original
						<INFIVE>; <INFIVE>; <INFIVE>; <INFIVE>; 	# skip next 4 lines in pb read file (equal to 1 read) -> erase this read
					}
					else{						# Read speichern
						$pbin=<INFIVE>; print OUTFOUR $pbin;		# write next 1/4 lines in pb read file (equal to 1 read)
						$pbin=<INFIVE>; print OUTFOUR $pbin;		# write next 2/4 lines in pb read file (equal to 1 read)
						$pbin=<INFIVE>; print OUTFOUR $pbin;		# write next 3/4 lines in pb read file (equal to 1 read)
						$pbin=<INFIVE>; print OUTFOUR $pbin;		# write next 4/4 lines in pb read file (equal to 1 read)
					}
				}
				<INFOUR>;
				<INFOUR>;
				<INFOUR>;
			}
		}

		open(OUTFIVE,">","Sim_Pipeline/$pre/$pre-filtered$r1") or die "Can't create Sim_Pipeline/$pre/$pre-filtered$r1";	# create file for filtered illumina r1 reads
		open(OUTSIX,">","Sim_Pipeline/$pre/$pre-filtered$r2") or die "Can't create Sim_Pipeline/$pre/$pre-filtered$r2";		# create file for filtered illumina r2 reads
		open(INSIX,"Sim_Pipeline/$pre/$pre$r1") or die "Can't open Sim_Pipeline/$pre/$pre$r1";				# input illumina r1 file with length
		open(INSEVEN,"Sim_Pipeline/$pre/$pre$r1") or die "Can't open Sim_Pipeline/$pre/$pre$r1";				# input illumina r1 reads for filtering
		open(INEIGHT,"Sim_Pipeline/$pre/$pre$r2") or die "Can't open Sim_Pipeline/$pre/$pre$r2";				# input illumina r2 reads for filtering

		while($r1line=<INSIX>){								# parse illumina r1 file
			chomp($r1line);
			if($r1line=~m/^@(.+)_(\d+)_\d+_.+_.+_.+/){			### wenn Referenz-Zeile ($1=ID, $2=read-start)
				if($2>$sizehash1{$1}){						### wenn read-start größer ist, als länge der sequenz mit entsprechender ID -> read löschen
					<INSEVEN>; <INSEVEN>; <INSEVEN>; <INSEVEN>;	# skip next 4 lines in illumina r1 read file (equal to 1 read)
					<INEIGHT>; <INEIGHT>; <INEIGHT>; <INEIGHT>;	# skip next 4 lines in illumina r2 read file (equal to 1 read)
				}
				else{									# Read speichern
					$r1in=<INSEVEN>; print OUTFIVE $r1in;		# write next 1/4 lines in illumina r1 read file (equal to 1 read)
					$r1in=<INSEVEN>; print OUTFIVE $r1in;		# write next 2/4 lines in illumina r1 read file (equal to 1 read)
					$r1in=<INSEVEN>; print OUTFIVE $r1in;		# write next 3/4 lines in illumina r1 read file (equal to 1 read)
					$r1in=<INSEVEN>; print OUTFIVE $r1in;		# write next 4/4 lines in illumina r1 read file (equal to 1 read)

					$r2in=<INEIGHT>; print OUTSIX $r2in;		# write next 1/4 lines in illumina r2 read file (equal to 1 read)
					$r2in=<INEIGHT>; print OUTSIX $r2in;		# write next 2/4 lines in illumina r2 read file (equal to 1 read)
					$r2in=<INEIGHT>; print OUTSIX $r2in;		# write next 3/4 lines in illumina r2 read file (equal to 1 read)
					$r2in=<INEIGHT>; print OUTSIX $r2in;		# write next 4/4 lines in illumina r2 read file (equal to 1 read)
				}
			
			}
		}
		
		close(OUTFOUR);
	}
	else{
		print OUTONE "$line1\n"		#write missing genome-files in file
	}
}

		


