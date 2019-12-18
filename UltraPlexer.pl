use strict;
use Data::Dumper;
use List::Util qw/max sum/;
use List::MoreUtils qw/mesh all/;
use Getopt::Long;   
use File::Path qw(make_path remove_tree);
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use Math::GSL::Randist qw/gsl_ran_binomial_pdf/;


##########################################################################
# MODIFY the following path to the Cortex directory - this script expects to
# find binaries compiled for 20, 40 and 60 colors.
my $cortex_bin_dir = '/gpfs/project/dilthey/software/CORTEX_release_v1.0.5.21/bin/';
##########################################################################

foreach my $colors (20, 40, 60)
{
	my $bin_fn = "${cortex_bin_dir}/cortex_var_31_c${colors}";
	die "Expected Cortex binary for $colors colors not found in directory $cortex_bin_dir - if necessary, modify this script so that the \$cortex_bin_dir variable points to right directory. Expected filename: $bin_fn" unless(-e $bin_fn);
}

$| = 1; 

# Debugging commands - ignore:
# die Dumper(observed_intersection_to_identity(16, 100, 99), 2, observed_intersection_to_identity_2(16, 100, 99));
# die observed_kmer_intersection_likelihood(16, 100, 0.90, 80);

my $action = 'classify';
my $samples_file = '';
my $cortex_temp_dir = $FindBin::Bin . '/cortex_temp';
my $k = 19;
my $oneSample_cortex_height = 20;
my $allSamples_cortex_height = 21;
my $p_seqError = 0.15;
my $prefix;
my $classificationSource_forCallFile = 'kmers';
my $testDataFromSimulation;
my $testDataFromRealBarcodedRun;
my $longReads_FASTQ;
my $ignoreAmbiguousReads;
GetOptions ( 
	'action:s' => \$action,
	'samples_file:s' => \$samples_file,
	'cortex_temp_dir:s' => \$cortex_temp_dir,
	'oneSample_cortex_height:s' => \$oneSample_cortex_height,
	'allSamples_cortex_height:s' => \$allSamples_cortex_height,
	'k:s' => \$k,
	'prefix:s' => \$prefix,
	'classificationSource:s' => \$classificationSource_forCallFile,
	'testDataFromSimulation:s' => \$testDataFromSimulation,
	'testDataFromRealBarcodedRun:s' => \$testDataFromRealBarcodedRun,
	'longReads_FASTQ:s' => \$longReads_FASTQ,
	'ignoreAmbiguousReads:s' => \$ignoreAmbiguousReads,
);

die "Please use --k <= 31" unless($k <= 31);

unless(-d $cortex_temp_dir)
{
	mkdir($cortex_temp_dir) or die "Cannot mkdir $cortex_temp_dir (you can specify a custom directory via --cortex_temp_dir)";
}

unless($prefix)
{
	die "Please specify --prefix (make sure this is unique per run)";
}	

if($action eq 'classifyTestData')
{
	unless($testDataFromSimulation or $testDataFromRealBarcodedRun)
	{
		die "You're in mode classifyTestData, but neither --testDataFromSimulation nor --testDataFromRealBarcodedRun is set to 1";
	}

	my $fn_output = $prefix . '.classification' . '_k' . $k;
	my $fn_output_OK = $prefix . '.classification' . '_k' . $k . '.done';
	if(-e $fn_output_OK)
	{
		die "File $fn_output_OK already existing - exit now, delete file if you want to redo.";
	}
	
	# die unless(-e $cortex_bin);^^
	die unless(-e $cortex_temp_dir);

	# build Cortex graphs and determine thresholds
	my @sampleIDs;
	my %longReads_per_sampleID;
	my %_ilmn_unique;
	my %basesPerSample;
	open(SAMPLES, '<', $samples_file) or die "Cannot open --samples_file $samples_file";
	while(<SAMPLES>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) == 5);
		my $sampleID = $line_fields[0];
		die "Sample IDs in $samples_file contain whitespaces -- can't deal with this. Offending example: $sampleID" if($sampleID =~ /\s/);
		my $nanoporeFile = $line_fields[1];
		die if($_ilmn_unique{$line_fields[2]});
		die if($_ilmn_unique{$line_fields[3]});
		my $Illumina_FASTQs_aref = [[$line_fields[2], $line_fields[3]]];
		$_ilmn_unique{$line_fields[2]} = 1;
		$_ilmn_unique{$line_fields[3]} = 1;
		$basesPerSample{$sampleID} = buildCortexFromFASTQ($prefix, $sampleID, $Illumina_FASTQs_aref, $k);
		push(@sampleIDs, $sampleID);
		$longReads_per_sampleID{$sampleID} = $nanoporeFile;
	}
	close(SAMPLES);
		
	my $readIDs_href = {};
	my $fh_output;
	open($fh_output, '>', $fn_output) or die "Cannot open $fn_output";
	print {$fh_output} join("\t", "readID", "readLength", "truthFromSimulation", "truthFromBarcodes", "maximumEstimatedIdentity", (map {'P' . $_ . '_adaptive'} 0 .. $#sampleIDs), (map {'P' . $_ . '_naive'} 0 .. $#sampleIDs), (map {'P' . $_ . '_proper'} 0 .. $#sampleIDs), (map {'kMers' . $_} 0 .. $#sampleIDs)), "\n";
	for(my $sampleI = 0; $sampleI <= $#sampleIDs; $sampleI++)
	{
		my $sampleID = $sampleIDs[$sampleI];
		classifyLongReads(
			$prefix,
			\@sampleIDs,
			$k,
			$longReads_per_sampleID{$sampleID},
			$fh_output,
			($testDataFromSimulation ? $sampleI : undef),
			($testDataFromRealBarcodedRun ? $sampleI : undef),
			$readIDs_href,
			$ignoreAmbiguousReads,
		);
	}
	close($fh_output);

	print "\nClassification done. Generated file $fn_output\n\n";
	
	open(OK, '>', $fn_output_OK) or die "Cannot open $fn_output_OK";
	print OK 1;
	close(OK);
}
elsif($action eq 'classify')
{
	unless((defined $longReads_FASTQ) and (-e $longReads_FASTQ))
	{
		die "Please specify a valid long-read FASTQ for classification via --longReads_FASTQ";
	}
	
	my $fn_output = $prefix . '.classification' . '_k' . $k;
	my $fn_output_OK = $prefix . '.classification' . '_k' . $k . '.done';
	if(-e $fn_output_OK)
	{
		die "File $fn_output_OK already existing - exit now, delete file if you want to redo.";
	}
	
	# die unless(-e $cortex_bin);
	mkdir($cortex_temp_dir) unless(-e $cortex_temp_dir);

	# build Cortex graphs and determine thresholds
	my @sampleIDs;
	my %longReads_per_sampleID;	my %_ilmn_unique;
	my %basesPerSample;
	open(SAMPLES, '<', $samples_file) or die "Cannot open --samples_file $samples_file";
	while(<SAMPLES>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die "Error - line $. in file $samples_file does not have three tab-delimited fields" unless(scalar(@line_fields) == 3);
		my $sampleID = $line_fields[0];
		die "Sample IDs in $samples_file contain whitespaces -- can't deal with this. Offending example: $sampleID" if($sampleID =~ /\s/);
		die if($_ilmn_unique{$line_fields[1]});
		die if($_ilmn_unique{$line_fields[2]});
		my $Illumina_FASTQs_aref = [[$line_fields[1], $line_fields[2]]];
		$_ilmn_unique{$line_fields[1]} = 1;
		$_ilmn_unique{$line_fields[2]} = 1;
		$basesPerSample{$sampleID} = buildCortexFromFASTQ($prefix, $sampleID, $Illumina_FASTQs_aref, $k);
		push(@sampleIDs, $sampleID);
	}
	close(SAMPLES);
	
	my $readIDs_href = {};
	my $fh_output;
	open($fh_output, '>', $fn_output) or die "Cannot open $fn_output";
	print {$fh_output} join("\t", "readID", "readLength", "truthFromSimulation", "truthFromBarcodes", "maximumEstimatedIdentity", (map {'P' . $_ . '_adaptive'} 0 .. $#sampleIDs), (map {'P' . $_ . '_naive'} 0 .. $#sampleIDs), (map {'P' . $_ . '_proper'} 0 .. $#sampleIDs), (map {'kMers' . $_} 0 .. $#sampleIDs)), "\n";
	classifyLongReads(
		$prefix,
		\@sampleIDs,
		$k,
		$longReads_FASTQ,
		$fh_output,
		undef,
		undef,
		$readIDs_href,
		$ignoreAmbiguousReads
	);
	close($fh_output);

	open(OK, '>', $fn_output_OK) or die "Cannot open $fn_output_OK";
	print OK 1;
	close(OK);
	
	print "\nClassification done. Generated file $fn_output\n\n";
} 
elsif($action eq 'generateCallFile')
{
	my $fn_input = $prefix . '.classification' . '_k' . $k;
	die "Classification file $fn_input not existing yet" unless(-e $fn_input);
	
	die "Please specify parameter --classificationSource, e.g. full, naive or kmers" unless(defined $classificationSource_forCallFile);
	my $fn_output = $fn_input . '.called_' . $classificationSource_forCallFile;
	
	my @sampleIDs;
	open(SAMPLES, '<', $samples_file) or die "Cannot open --samples_file $samples_file";
	while(<SAMPLES>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless((scalar(@line_fields) == 5) || (scalar(@line_fields) == 3));
		my $sampleID = $line_fields[0];
		push(@sampleIDs, $sampleID);
	}
	close(SAMPLES);

	open(OUT, '>', $fn_output) or die "Cannot open $fn_output";
	print "Reading $fn_input in classification mode $classificationSource_forCallFile...\n";
	open(INPUT, '<', $fn_input) or die "Cannot open $fn_input";
	my $headerLine = <INPUT>;
	chomp($headerLine);
	my @header_fields = split(/\t/, $headerLine);
	my %readLengths;
	while(<INPUT>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die "Field number mismatch in line $. of $fn_input -- $#line_fields vs $#header_fields" unless($#line_fields == $#header_fields);
		my %line = (mesh @header_fields, @line_fields);
		my $readID = $line{readID}; die unless(defined $readID);
		my $readLength = $line{readLength}; die unless(defined $readLength);
		
		die "Read length mismatch" unless((not defined $readLengths{$readID}) or ($readLengths{$readID} == $readLength));
		$readLengths{$readID} = $readLength;
		
		my @callDistribution = getCallDistributionFromLine(\%line, \@sampleIDs, $classificationSource_forCallFile);
		
		my $maximumEstimatedIdentity = $line{maximumEstimatedIdentity};
		die unless(defined $maximumEstimatedIdentity);
		
		my $call = callDistribution(\@callDistribution, $ignoreAmbiguousReads);
		my $classification = $sampleIDs[$call];
			
		print OUT join("\t", $readID, $classification, $maximumEstimatedIdentity), "\n";
	}
		
	close(OUT);
	
	print "\nGenerated call file: $fn_output \n\n";
}
elsif($action eq 'evaluate')
{
	my $fn_input = $prefix . '.classification' . '_k' . $k;
	die "Classification file $fn_input not existing yet" unless(-e $fn_input);
	
	my @sampleIDs;
	open(SAMPLES, '<', $samples_file) or die "Cannot open --samples_file $samples_file";
	while(<SAMPLES>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) == 5);
		my $sampleID = $line_fields[0];
		push(@sampleIDs, $sampleID);
	}
	close(SAMPLES);
	my %sampleID_2_i;
	for(my $i = 0; $i <= $#sampleIDs; $i++)
	{
		$sampleID_2_i{$sampleIDs[$i]} = $i;
	}
	
	(my $read_mapping_truth_href, my $read_mapping_distances_abs_href, my $read_mapping_distances_prop_href) = readReadTruthFromMapping(); # todo 
	my $read_simulationOrBarcodes_href = readReadTruthFromSimulationOrBarcodes($fn_input, \@sampleIDs); # todo
	#my $read_truth_href = {};
	#my $read_barcodes_href = {};
	
	my %readClassification_distribution;
	my %readLengths;
	my %maximumIdentities;
	#foreach my $classificationSource (qw/adaptive simulationOrBarcodes full naive kmers/) # todo
	my @evaluateMethods = qw/random full naive kmers adaptive/;
	foreach my $classificationSource (@evaluateMethods) 
	{
		print "Reading $fn_input in classification mode $classificationSource...\n";
		open(INPUT, '<', $fn_input) or die "Cannot open $fn_input";
		my $headerLine = <INPUT>;
		chomp($headerLine);
		my @header_fields = split(/\t/, $headerLine);
		while(<INPUT>)
		{
			# last if($. > 2000); # todos
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line, -1);
			die "Field number mismatch in line $. of $fn_input -- $#line_fields vs $#header_fields" unless($#line_fields == $#header_fields);
			my %line = (mesh @header_fields, @line_fields);
			my $readID = $line{readID}; die unless(defined $readID);
			my $readLength = $line{readLength}; die unless(defined $readLength);
			my $truth = extractSimulationOrBarcodTruth(\%line); die unless(defined $truth);
			
			die Dumper("Read length mismatch", $readID, $readLength, $readLengths{$readID}) unless((not defined $readLengths{$readID}) or ($readLengths{$readID} == $readLength));
			$readLengths{$readID} = $readLength;
			
			my $maximumEstimatedIdentity = $line{maximumEstimatedIdentity};
			die unless(defined $maximumEstimatedIdentity);	
			die Dumper("Maximum identity mismatch", $readID, $maximumEstimatedIdentity, $maximumIdentities{$readID}) unless((not defined $maximumIdentities{$readID}) or ($maximumIdentities{$readID} == $maximumEstimatedIdentity));
			$maximumIdentities{$readID} = $maximumEstimatedIdentity;
			
			if($classificationSource eq 'random')
			{
				my @callDistribution = getRandomCallDistribution(\@sampleIDs);		
				$readClassification_distribution{$classificationSource}{$readID} = \@callDistribution;
			}
			elsif($classificationSource ne 'simulationOrBarcodes')
			{
				my @callDistribution = getCallDistributionFromLine(\%line, \@sampleIDs, $classificationSource);
				$readClassification_distribution{$classificationSource}{$readID} = \@callDistribution;
			}
			else
			{
				if($truth <= $#sampleIDs)
				{
					my @callDistribution_barcodes = ((0) x scalar(@sampleIDs));
					die unless(scalar(@callDistribution_barcodes) == scalar(@sampleIDs));
					$callDistribution_barcodes[$truth] = 1;
					$readClassification_distribution{'simulationOrBarcodes'}{$readID} = \@callDistribution_barcodes;	
				}
				else
				{
					die Dumper("Something weird has happened here if this is simulated data!"); 
				}
			}
		}
		close(INPUT);
	}

	my $fn_calibration_output = $fn_input . '.calibrationAgainstTruthOrBarcodes';
	open(CALIBRATION, '>', $fn_calibration_output) or die;
	foreach my $classificationSource (sort keys %readClassification_distribution)
	{
		print "Evaluatin calibration $classificationSource...\n";
		my $callis_1 = 0;
		foreach my $readID (keys %{$readClassification_distribution{$classificationSource}})
		{
			foreach my $sampleID (@sampleIDs)
			{
				my $sampleID_i = $sampleID_2_i{$sampleID};
				die unless(defined $sampleID_i);
				my $p_sampleID = $readClassification_distribution{$classificationSource}{$readID}[$sampleID_i];
				my $isCorrect = ($read_simulationOrBarcodes_href->{$readID}{$sampleID}) ? 1 : 0;
				print CALIBRATION join("\t", $classificationSource, $readID, $sampleID, $p_sampleID, $isCorrect), "\n";
			}
		}
	}
	
	close(CALIBRATION);
	
	
	my %readClassification_calls;	
	foreach my $classificationSource (sort keys %readClassification_distribution)
	{
		print "Calling $classificationSource...\n";
		my $callis_1 = 0;
		foreach my $readID (keys %{$readClassification_distribution{$classificationSource}})
		{
			my $call = callDistribution($readClassification_distribution{$classificationSource}{$readID}, $ignoreAmbiguousReads);
			# print Dumper($readClassification_distribution{$classificationSource}{$readID}, $call); # todo
			$readClassification_calls{$classificationSource}{$readID} = $sampleIDs[$call];
			# if($call == 0)
			# {
				# $callis_1++;
			# }
		}
		# print "$classificationSource callis_1: $callis_1\n";			
	}
	
	# exit; # todo
	
	# evaluation against proper one-dimensional truth
	
	print "Accuracy by truth sample:\n";	
	foreach my $classificationSource (@evaluateMethods)
	{
		my %reads_correct_by_truth;
		my %read_distances_abs_byTruth;
		my %read_distances_prop_byTruth;

		foreach my $readID (keys %{$readClassification_calls{$classificationSource}})
		{
			die unless(exists $read_mapping_truth_href->{$readID});
			die unless(exists $read_simulationOrBarcodes_href->{$readID});
			my @truth_samples = keys %{$read_simulationOrBarcodes_href->{$readID}};
			die unless(scalar(@truth_samples) == 1);
			my $truth_sample = $truth_samples[0];
			my $assignment = $readClassification_calls{$classificationSource}{$readID}; die unless(defined $assignment);
			$reads_correct_by_truth{$truth_sample}[0]++;
			if($truth_sample eq $assignment)
			{
				$reads_correct_by_truth{$truth_sample}[1]++;			
			}
			push(@{$read_distances_abs_byTruth{$truth_sample}}, $read_mapping_distances_abs_href->{$readID}{$assignment});
			push(@{$read_distances_prop_byTruth{$truth_sample}}, $read_mapping_distances_prop_href->{$readID}{$assignment});
		}
		
		print "\tSource $classificationSource\n";
		for(my $sampleI = 0; $sampleI <= $#sampleIDs; $sampleI++)
		{
			my $truth = $sampleIDs[$sampleI];
			die unless(scalar(@{$read_distances_abs_byTruth{$truth}}));
			die unless(scalar(@{$read_distances_prop_byTruth{$truth}}));
			print "\t\t$truth\n";
			print "\t\t\tAbsolute classifcation accuracy: $reads_correct_by_truth{$truth}[0] $reads_correct_by_truth{$truth}[1] ", sprintf("%.2f", $reads_correct_by_truth{$truth}[1] / $reads_correct_by_truth{$truth}[0]), "\n"; 
			print "\t\t\tAbsolute distance to target: ", getMean($read_distances_abs_byTruth{$truth}, 3), " (mean) - ", getMedian($read_distances_abs_byTruth{$truth}, 3), " (median).\n";
			print "\t\t\tRelative distance to target: ", getMean($read_distances_prop_byTruth{$truth}, 3), " (mean) - ", getMedian($read_distances_prop_byTruth{$truth}, 3), " (median).\n";
			
			if($truth eq 'GCF_000013465.1_ASM1346v1_genomic')
			{
				# die Dumper([sort @{$read_distances_abs_byTruth{$truth}}]);
			}	
		}

		
	}
	# foreach my $evaluationData ([$fn_input . '.evaluationAgainstMapping', $read_mapping_truth_href], [$fn_input . '.evaluationAgainstSimulationOrBarcodes', $read_simulationOrBarcodes_href])
	foreach my $evaluationData ([$fn_input . '.evaluationAgainstMapping', $read_mapping_truth_href])
	{
		my $fn_evaluation_output = $evaluationData->[0];
		my $truthSet = $evaluationData->[1];

		my $fh_evaluation_output;
		open($fh_evaluation_output, '>', $fn_evaluation_output) or die "Cannot open $fn_evaluation_output";
		print {$fh_evaluation_output} join("\t", qw/readID readLength maximumIdentity method correct classifiedAs truth distanceToTruthAbs distanceToTruthRel/), "\n";
		foreach my $classificationSource (@evaluateMethods)
		{	
			evaluateReadAssignments($readClassification_calls{$classificationSource}, $truthSet, $read_mapping_distances_abs_href, $read_mapping_distances_prop_href, \%readLengths, \%maximumIdentities, $classificationSource, $fh_evaluation_output);
		}
		close($fh_evaluation_output);
		print "\n\nGenerated file: $fn_evaluation_output\n\n";
	}
}


sub evaluateReadAssignments
{
	my $assignments_href = shift;
	my $truth_href = shift;
	my $read_mapping_distances_abs_href = shift;
	my $read_mapping_distances_prop_href = shift;
	
	my $readLengths_href = shift;
	my $maximumIdentities_href = shift;
	my $label = shift;
	my $fh_evaluation_output = shift;
	
	my $total_reads = 0;
	my $total_noTruth = 0;
	my $total_reads_OK = 0;
	my %reads_per_source;
	my %reads_per_source_correct;
	my %reads_false_calls_to;
	my %truth_reads;
	foreach my $readID (keys %$assignments_href)
	{
		unless(exists $truth_href->{$readID})
		{
			#warn "No truth for $readID";
			$total_noTruth++;
			next;
		}
		$reads_per_source{$assignments_href->{$readID}}++;
		
		my $isCorrect = 0;
		$total_reads++;
		
		if($truth_href->{$readID}{$assignments_href->{$readID}})
		{	
			$total_reads_OK++;
			$isCorrect = 1;
		}
		else
		{
			if($label eq 'kmers')
			{
				# print join("\t", $readID, $assignments_href->{$readID}, join(';', keys %{$truth_href->{$readID}})), "\n";
			}	
			
			$reads_false_calls_to{$assignments_href->{$readID}}++;
		}
		
		foreach my $key_truth (keys %{$truth_href->{$readID}})
		{
			$truth_reads{$key_truth}++;
		}
		
		die unless(defined $readLengths_href->{$readID});
		print {$fh_evaluation_output} join("\t", $readID, $readLengths_href->{$readID}, $maximumIdentities_href->{$readID}, $label, $isCorrect, $assignments_href->{$readID}, join(';', keys %{$truth_href->{$readID}}), $read_mapping_distances_abs_href->{$readID}{$assignments_href->{$readID}}, $read_mapping_distances_prop_href->{$readID}{$assignments_href->{$readID}}), "\n";
		
	}
	
	print "Read assignment evaluation for $label\n";
	print "\t", "Total reads: ", $total_reads, "\n";
	print "\t", "Assigned OK: ", $total_reads_OK, " (", sprintf("%.2f", $total_reads_OK/$total_reads * 100), "%)\n";
	print "\t", "No truth   : ", $total_noTruth, "\n";
	print "\n";
	print "\tReads come from (assigned truth; multiples are counted!):\n";
	foreach my $source (sort keys %truth_reads)
	{
		print "\t\t", $source, ": ", $truth_reads{$source}, " reads.\n";
	}	
	
	print "\tAssigned to:\n";
	foreach my $source (sort keys %reads_per_source)
	{
		print "\t\t", $source, ": ", $reads_per_source{$source}, " reads.\n";
	}
	print "\tFalse calls go to:\n";
	foreach my $target (sort keys %reads_false_calls_to)
	{
		print "\t\t", $target, ": ", $reads_false_calls_to{$target}, " reads.\n";
	}	
	print "\n";
}

sub getMean
{
	my $aref = shift;
	my $digits = shift;
	my $S = sum(@$aref);
	die "getMean(..): don't call me with empty input - " . scalar(@$aref) unless(scalar(@$aref));
	my $mean = $S / scalar(@$aref);
	$mean = sprintf("%." . $digits . "f", $mean) if($digits);
	return $mean;
}

sub getMedian
{
	my $aref = shift;
	my $digits = shift;
	my @a_sorted = sort(@$aref);
	die unless(@a_sorted);
	my $median = $a_sorted[$#a_sorted / 2];
	$median = sprintf("%." . $digits . "f", $median) if($digits);
	return $median;
}

sub readReadTruthFromMapping
{
	my $fn = $samples_file . '.readStats';

	my %forReturn;
	my %forReturn_distances_abs;
	my %forReturn_distances_prop;
	
	open(F, '<', $fn) or die "Cannot open $fn";
	my $headerLine = <F>;
	chomp($headerLine);
	my @header_fields = split(/\t/, $headerLine);
	my @sampleIDs_in_truth = @header_fields[3 .. $#header_fields];
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless($#line_fields == $#header_fields);
		my %line = (mesh @header_fields, @line_fields);
		my $readID = $line{readID}; die unless(defined $readID);
		
		my $max_aligned_bases;
		foreach my $sampleID (@sampleIDs_in_truth)
		{
			die unless(defined $line{$sampleID});
			if((not defined $max_aligned_bases) or ($max_aligned_bases < $line{$sampleID}))
			{
				$max_aligned_bases = $line{$sampleID};
			}
		}
		die unless(defined $max_aligned_bases);
		my $readLength = $line{readLength};
		die unless($readLength);	
		# next unless($readLength == $max_aligned_bases);
		die "File $fn line $., rL $readLength, max alignment $max_aligned_bases" unless($readLength == $max_aligned_bases);

		foreach my $sampleID (@sampleIDs_in_truth)
		{
			die unless(defined $line{$sampleID});
			if($line{$sampleID} == $max_aligned_bases)
			{
				$forReturn{$readID}{$sampleID}++;
			}

			my $distance_abs = $max_aligned_bases - $line{$sampleID};
			my $distance_prop = $distance_abs / $readLength;

			$forReturn_distances_abs{$readID}{$sampleID} = $distance_abs;
			$forReturn_distances_prop{$readID}{$sampleID} = $distance_prop;
		}	

		die unless(defined $forReturn{$readID});	
	}
	close(F);
	
	return (\%forReturn, \%forReturn_distances_abs, \%forReturn_distances_prop);
}

	
sub getCallDistributionFromLine
{
	my $line_href = shift;
	my $sampleIDs_aref = shift;
	my $classificationSource = shift;
	
	my @callDistribution;
	for(my $i = 0; $i <= $#{$sampleIDs_aref}; $i++)
	{
		my $v;
		if($classificationSource eq 'full')
		{
			$v = $line_href->{'P' . $i . '_proper'};
		}
		elsif($classificationSource eq 'naive')
		{
			$v = $line_href->{'P' . $i . '_naive'};				
		}
		elsif($classificationSource eq 'adaptive')
		{
			$v = $line_href->{'P' . $i . '_adaptive'};				
		}		
		elsif($classificationSource eq 'kmers')
		{
			$v = $line_href->{'kMers' . $i};								
		}
		elsif($classificationSource eq 'random')
		{
			$v = 1/scalar(@{$sampleIDs_aref});
		}
		else
		{
			die "Unknown classification source: $classificationSource";
		}
		die unless(defined $v);
		push(@callDistribution, $v);
	}
	
	if($classificationSource eq 'kmers')
	{
		@callDistribution = kMerVectorToDistribution(@callDistribution);
	}
	
	return @callDistribution;
				
}

sub getRandomCallDistribution
{
	my $sampleIDs_aref = shift;
	my $selection = int(rand(scalar(@$sampleIDs_aref)));
	
	my @callDistribution = ((0) x scalar(@$sampleIDs_aref));
	die unless($#{$sampleIDs_aref} == $#callDistribution);
	
	die unless($selection >= 0);
	die unless($selection <= $#callDistribution);
	$callDistribution[$selection] = 1;
	
	return @callDistribution;
}

sub extractSimulationOrBarcodTruth
{
	my $line_href = shift;
	my $truthFromSimulation = $line_href->{truthFromSimulation};
	my $truthFromBarcodes = $line_href->{truthFromBarcodes};
	$truthFromSimulation = undef if($truthFromSimulation eq 'NA');
	$truthFromBarcodes = undef if($truthFromBarcodes eq 'NA');
	
	die Dumper("Weird line", $line_href) unless((defined $truthFromSimulation) xor (defined $truthFromBarcodes));
	
	my $truth = ((defined $truthFromSimulation) ? $truthFromSimulation : $truthFromBarcodes);	
	
	return $truth;
}

sub readReadTruthFromSimulationOrBarcodes
{
	my $fn_input = shift;
	my $sampleIDs_aref = shift;
	
	my %forReturn;
	
	open(INPUT, '<', $fn_input) or die "Cannot open $fn_input";
	my $headerLine = <INPUT>;
	chomp($headerLine);
	my @header_fields = split(/\t/, $headerLine);
	while(<INPUT>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die "Field number mismatch in line $. of $fn_input -- $#line_fields vs $#header_fields" unless($#line_fields == $#header_fields);
		my %line = (mesh @header_fields, @line_fields);
		my $readID = $line{readID}; die unless(defined $readID);
		my $readLength = $line{readLength}; die unless(defined $readLength);
		my $truth = extractSimulationOrBarcodTruth(\%line);
		my $sampleID_truth = $sampleIDs_aref->[$truth];
		die unless(defined $sampleID_truth);
		
		$forReturn{$readID}{$sampleID_truth} = 1;
	}
	close(INPUT);
	
	return \%forReturn;
}

sub callDistribution
{
	my $distribution_aref = shift;
	my $ignoreAmbiguousReads = shift;
	
	my $distr_cum = 0;
	my $distr_max;
	for(my $i = 0; $i <= $#{$distribution_aref}; $i++)
	{
		$distr_cum += $distribution_aref->[$i];
		if((not defined $distr_max) or ($distr_max < $distribution_aref->[$i]))
		{
			$distr_max = $distribution_aref->[$i];
		}
	}
	die unless(defined $distr_max);
	my @distr_elements_with_max = grep {$_ == $distr_max} @$distribution_aref;
	
	if($ignoreAmbiguousReads and (scalar(@distr_elements_with_max) > 1))
	{
		return 'NA';
	}
	die Dumper("Invalid distribution", $distribution_aref, $distr_cum) unless(abs(1 - $distr_cum) <= 1e-6);
	my $v = rand(1);
	
	$distr_cum = 0;
	for(my $i = 0; $i <= $#{$distribution_aref}; $i++)
	{
		$distr_cum += $distribution_aref->[$i];
		$distr_cum = 1 if($i == $#{$distribution_aref});
		if($v <= $distr_cum)
		{
			return $i;
		}
	}
	
	die "This should not happen";
}

sub kMerVectorToDistribution
{
	my @input = @_;
	my $max;
	my $max_which;
	for(my $i = 0; $i <= $#input; $i++)
	{
		if((not defined $max) or ($max < $input[$i]))
		{
			$max = $input[$i];
			$max_which = $i;
		}
	}
	
	my $max_howMany = 0;
	for(my $i = 0; $i <= $#input; $i++)
	{
		if($input[$i] == $max)
		{
			$max_howMany++;
		}
	}	
	
	for(my $i = 0; $i <= $#input; $i++)
	{
		if($input[$i] == $max)
		{
			$input[$i] = 1/$max_howMany;
		}
		else
		{
			$input[$i] = 0;
		}
	}		
	
	return @input;
}

sub log_likelihoods_to_posterior_probabilities
{
	my @ll = @_;

	my $max_ll = max(@ll);
	foreach my $ll (@ll)
	{
		$ll -= $max_ll;
	}
	
	
	foreach my $ll (@ll)
	{
		$ll = exp($ll);
	}
	

	die unless(all {($_ >= 0) and ($_ <= 1)} @ll);
	die unless((max @ll) > 0);
			
	my $s_l = sum(@ll);
		
	foreach my $ll (@ll)
	{
		$ll = $ll/$s_l;
	}
		
	die unless(abs(1 - sum(@ll)) <= 1e-3);
	
	return @ll;
}

sub classifyLongReads
{
	my $prefix = shift;
	my $sampleIDs_aref = shift;
	my $k = shift;
	my $FASTQ_longReads = shift;
	my $fh_output = shift;
	my $truthColour_fromSimulation = shift;
	my $truthColour_fromBarcodes = shift;
	my $seenReadIDs_href = shift;
	my $classifyLongReads = shift;
	die unless(defined $seenReadIDs_href);
	die unless(defined $classifyLongReads);
	
	$truthColour_fromSimulation = 'NA' unless(defined $truthColour_fromSimulation);
	$truthColour_fromBarcodes = 'NA' unless(defined $truthColour_fromBarcodes);
	
	die unless(-e $FASTQ_longReads);
	
	my $longestRead_length = getLongestReadLength_FASTQ($FASTQ_longReads);
	my $n_reads_logReadsFASTQ = countReadsInFASTQ($FASTQ_longReads);
	
	my $kMer_survival_rate = (1-$p_seqError)**$k;
	my $n_total_effective_kmers = (4**$k)/2;
	
	my @n_kmers_above_threshold;
	my @coverage_of_kmers_above_threshold;
	my @kmer_threshold;
	my @logP_novelkMers;
	my %sampleID_2_colour;
	
	my $fn_output_meta = $cortex_temp_dir . '/meta_' . $prefix . '_' . $k;	
	open(META, '>', $fn_output_meta) or die "Cannot open $fn_output_meta";	
	print META join("\t", qw/sampleI sampleID kMerThreshold kMersBelowT kMersAboveT coverageOfkMersBelowT coverageOfkMersAboveT logPnovelkMer/), "\n";

	my $fn_binaries_list = $cortex_temp_dir . '/binariesList_' . $prefix . '_' . $k;
	open(FOFN_OUTER, '>', $fn_binaries_list) or die "Cannot open $fn_binaries_list";
	for(my $colourI = 0; $colourI <= $#{$sampleIDs_aref}; $colourI++)
	{
		my $fn_binaries_list_thisColour = $cortex_temp_dir . '/binariesList_' . $prefix . '_' . $k . '_c' . $colourI;

		my $sampleID = $sampleIDs_aref->[$colourI];
		my $cortexBinary = getCortexPathForSampleID($prefix, $sampleID, $k, 1);
		
		open(FOFN_INNER, '>', $fn_binaries_list_thisColour) or die "Cannot open $fn_binaries_list_thisColour";		
		print FOFN_INNER abs_path($cortexBinary), "\n";
		close(FOFN_INNER);
		
		print FOFN_OUTER abs_path($fn_binaries_list_thisColour), "\n";
		
		$sampleID_2_colour{$sampleID} = $colourI;		
		
		(my $n_kmers_belowT, my $n_kmers_aboveT) = getkMerThresholdStats($prefix, $sampleID, $k);
		(my $coverage_kmers_belowT, my $coverage_kmers_aboveT) = getkMerThresholdCoverageStats($prefix, $sampleID, $k);

		push(@n_kmers_above_threshold, $n_kmers_aboveT);
		push(@coverage_of_kmers_above_threshold, $coverage_kmers_aboveT);
		push(@kmer_threshold, getkMerThreshold($prefix, $sampleID, $k));
		
		my $logP_novelkMer = log(10 * (1/($n_total_effective_kmers - $n_kmers_above_threshold[$colourI])));
		my $logP_presentkMer_minimum = log($kmer_threshold[$colourI]/$coverage_of_kmers_above_threshold[$colourI]);
		warn "Sample $sampleID (prefix $prefix), logP_novelkMer = $logP_novelkMer, and logP_presentkMer_minimum = $logP_presentkMer_minimum" unless($logP_presentkMer_minimum > $logP_novelkMer);
		push(@logP_novelkMers, $logP_novelkMer);
		
		print META join("\t",
			$colourI,
			$sampleID,
			$kmer_threshold[$#kmer_threshold],
			$n_kmers_belowT,
			$n_kmers_aboveT,
			$coverage_kmers_belowT,
			$coverage_kmers_aboveT,
			$logP_novelkMers[$#logP_novelkMers],
		), "\n";
	}
	close(FOFN_OUTER);
	close(META);
	
	my $FASTA_longReads = $FASTQ_longReads . '.' . $prefix . '.fasta';
	fastq_to_fasta($FASTQ_longReads, $FASTA_longReads);
	
	my $fn_reads_list = $cortex_temp_dir . '/readsList_' . $prefix . '_' . $k;
	open(RLIST, '>', $fn_reads_list) or die "Cannot open $fn_reads_list";
	print RLIST abs_path($FASTA_longReads), "\n";
	close(RLIST);
	
	my $capture_output_file = $FASTA_longReads . '.colour_covgs.output';

	my $use_max_read_len = 1.5 * $longestRead_length;
	my $cortex_bin = get_cortex_bin_by_nColours(scalar(@$sampleIDs_aref));
	my $cortex_cmd = qq(/usr/bin/time -v $cortex_bin --mem_height $allSamples_cortex_height --mem_width 100 --colour_list $fn_binaries_list --kmer_size ${k} --align ${fn_reads_list},no --align_input_format LIST_OF_FASTA --max_read_len $use_max_read_len &> $capture_output_file);
	system($cortex_cmd) and die "Command $cortex_cmd failed -- check $capture_output_file for error messages";
	
	my $expected_coverages_file = $FASTA_longReads . '.colour_covgs';
	die "File $expected_coverages_file not existing" unless(-e $expected_coverages_file);
	
	my $currentReadID;
	my $currentReadLength;
	my @currentReadLines;
	
	my %thisFile_reads;
	my $processReadLines = sub {
		my $readID = shift;
		my $readLength = shift;
		
		if($readID)
		{
			die "I've already come across read ID '$readID' (length $readLength) (now processing $FASTQ_longReads) - please make sure that all input long read IDs are unique!"	if($seenReadIDs_href->{$readID});
			$seenReadIDs_href->{$readID}++;
			$thisFile_reads{$readID}++;
			
			my $debug = 0;
			# the following conditional catches sporadic mis-formatted lines (produced by Cortex)
			if(scalar(@currentReadLines) == scalar(@$sampleIDs_aref))
			{
				print "Classifying read $readID\n" if($debug);

				my $n_kmers;
				my @log_likelihoods_naive;
				my @log_likelihoods_proper;
				my @kmers_per_colour;
				
				for(my $colourI = 0; $colourI <= $#{$sampleIDs_aref}; $colourI++)
				{
					my $n_kMers_inSample = 0;
					my $n_kMers_notInSample = 0;
					
					my $ll_kMer_notPresent_naive = log(1e-8);
					my $ll_kMer_notPresent_proper = log(1-$kMer_survival_rate) + $logP_novelkMers[$colourI];
					
					my $ll_naive = 0;
					my $ll_proper = 0;
					my @kmer_coverages = split(/ /, $currentReadLines[$colourI]);
					
					if(defined $n_kmers)
					{
						die unless ($n_kmers == scalar(@kmer_coverages)); 
					}
					$n_kmers = scalar(@kmer_coverages) unless(defined $n_kmers);
					
					# warn Dumper($colourI, ['not present', log(1e-8)], ['present', log(1/$n_kmers_above_threshold[$colourI])]);
					foreach my $kMerCoverage (@kmer_coverages)
					{
						die unless($kMerCoverage =~ /^\d+$/);

						if($kMerCoverage >= $kmer_threshold[$colourI])
						{
							$ll_naive += log(1/$n_kmers_above_threshold[$colourI]);
							$ll_proper += (log($kMer_survival_rate) +  log($kMerCoverage/$coverage_of_kmers_above_threshold[$colourI]));
							$n_kMers_inSample++;
						}
						else
						{
							$ll_naive += $ll_kMer_notPresent_naive;
							$ll_proper += $ll_kMer_notPresent_proper;
							$n_kMers_notInSample++;
						}
					}
					push(@log_likelihoods_naive, $ll_naive);
					push(@log_likelihoods_proper, $ll_proper);
					push(@kmers_per_colour, $n_kMers_inSample);
					
					print "\tColour $colourI: $n_kMers_inSample/$n_kmers present, $n_kMers_notInSample/$n_kmers not present\n" if($debug);
				}
				die unless($n_kmers == ($readLength - $k + 1));
				
				my @log_likelihoods_adaptive;
				my $max_kMers_shared = max(@kmers_per_colour);
				my $maximumIdentity = observed_intersection_to_identity($k, $n_kmers, $max_kMers_shared);
				for(my $colourI = 0; $colourI <= $#{$sampleIDs_aref}; $colourI++)
				{
					push(@log_likelihoods_adaptive, observed_kmer_intersection_likelihood($k, $n_kmers, $maximumIdentity, $kmers_per_colour[$colourI]));
				}
				die unless(max(@log_likelihoods_adaptive) ==  observed_kmer_intersection_likelihood($k, $n_kmers, $maximumIdentity, $max_kMers_shared));
				
				my @posterior_probabilities_adaptive = log_likelihoods_to_posterior_probabilities(@log_likelihoods_adaptive);
				my @posterior_probabilites_naive = log_likelihoods_to_posterior_probabilities(@log_likelihoods_naive);
				my @posterior_probabilites_proper = log_likelihoods_to_posterior_probabilities(@log_likelihoods_proper);
				
				print {$fh_output} join("\t", $readID, $readLength, $truthColour_fromSimulation, $truthColour_fromBarcodes, $maximumIdentity, @posterior_probabilities_adaptive, @posterior_probabilites_naive, @posterior_probabilites_proper, @kmers_per_colour), "\n";
			}
		}
	};
	
	open(COV, '<', $expected_coverages_file) or die "Cannot open $expected_coverages_file";	
	while(<COV>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die "Weird line (I) $. in file $expected_coverages_file: $line" unless(substr($line, 0, 1) eq '>');
		
		my $thisLine_readID = substr($line, 1);
		$thisLine_readID =~ s/_colour_\d+_kmer_coverages//;
		
		if($thisLine_readID ne $currentReadID)
		{
			$processReadLines->($currentReadID, $currentReadLength);
			@currentReadLines = ();
			$currentReadID = $thisLine_readID;
			die "Weird line (II) $. in file $expected_coverages_file: $line ($currentReadID $thisLine_readID)" if($line =~ /kmer_coverages/);
			my $read_sequence = <COV>;
			chomp($read_sequence);
			$currentReadLength = length($read_sequence);
			die unless($currentReadLength <= $longestRead_length);
		}
		else
		{
			die "Weird unexpected line (III) $. in $expected_coverages_file: $line ($currentReadID $thisLine_readID)" unless($line =~ /_colour_(\d+)_kmer_coverages/);
			my $colour = $1;
			my $coverages = <COV>;
			chomp($coverages);
			if($colour <= $#{$sampleIDs_aref})
			{
				die unless($colour == scalar(@currentReadLines));
				
				push(@currentReadLines, $coverages);
			}
		}
	}
	close(COV);
	if($currentReadID)
	{
		$processReadLines->($currentReadID, $currentReadLength);
	}
	
	unlink($fn_binaries_list);	
	unlink($fn_reads_list);	
	unlink($FASTA_longReads); 
	unlink($expected_coverages_file);
	
	print "File $FASTQ_longReads: Classified ", scalar(keys %thisFile_reads), " of $n_reads_logReadsFASTQ long reads\n";
}

sub getkMerThresholdStats
{
	my $prefix = shift;
	my $sampleID = shift;
	my $k = shift;	
	die unless(defined $k);
	
	my $cortexFile = getCortexPathForSampleID($prefix, $sampleID, $k);
	my $coverage = $cortexFile . '.coverage';

	my $T = getkMerThreshold($prefix, $sampleID, $k);
	
	my $kMers_belowT = 0;
	my $kMers_aboveT = 0;
	open(C, '<', $coverage) or die "Cannot open $coverage";
	while(<C>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		die unless(scalar(@f) == 2);
		if($. == 0)
		{
			die unless($f[0] eq 'KMER_COVG');
			die unless($f[1] eq 'FREQUENCY');
		}
		else
		{
			$kMers_belowT += $f[1] if($f[0] < $T);
			$kMers_aboveT += $f[1] if($f[0] >= $T);
			
		}
	}
	close(C);
	
	print "Sample $sampleID: $kMers_belowT $k-mers below threshold of $T, $kMers_aboveT above (which will be used)\n";
	
	return ($kMers_belowT, $kMers_aboveT);
}

sub getkMerThresholdCoverageStats
{
	my $prefix = shift;
	my $sampleID = shift;
	my $k = shift;	
	die unless(defined $k);
	
	my $cortexFile = getCortexPathForSampleID($prefix, $sampleID, $k);
	my $coverage = $cortexFile . '.coverage';

	my $T = getkMerThreshold($prefix, $sampleID, $k);
	
	my $kMers_belowT = 0;
	my $kMers_aboveT = 0;
	my $kMers_coverage_belowT = 0;
	my $kMers_coverage_aboveT = 0;
	
	open(C, '<', $coverage) or die "Cannot open $coverage";
	while(<C>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		die unless(scalar(@f) == 2);
		if($. == 0)
		{
			die unless($f[0] eq 'KMER_COVG');
			die unless($f[1] eq 'FREQUENCY');
		}
		else
		{
			$kMers_belowT += $f[1] if($f[0] < $T);
			$kMers_aboveT += $f[1] if($f[0] >= $T);
			
			$kMers_coverage_belowT += ($f[0] * $f[1]) if($f[0] < $T);
			$kMers_coverage_aboveT += ($f[0] * $f[1]) if($f[0] >= $T);			
		}
	}
	close(C);
	
	print "Sample $sampleID: $kMers_belowT $k-mers below threshold of $T; coverage $kMers_coverage_belowT, $kMers_aboveT above (which will be used); coverage $kMers_coverage_aboveT\n";
	
	return ($kMers_coverage_belowT, $kMers_coverage_aboveT);
}


sub buildCortexFromFASTQ
{
	my $prefix = shift;
	my $sampleID = shift;
	my $FASTQs_aref = shift;
	my $k = shift;
	die unless(defined $k);
	
	my $totalBases = 0;
	foreach my $FASTQpair (@$FASTQs_aref)
	{
		die unless(scalar(@$FASTQpair) == 2);
		foreach my $FASTQ (@$FASTQpair)
		{
			die if($FASTQpair->[0] eq $FASTQpair->[1]);
			$totalBases += countBasesInFASTQ($FASTQ);		
		}
	}
	
	my $totalBases_Mb = sprintf("%.2f", $totalBases / (1000**2));
	print "Sample $sampleID, have $totalBases_Mb Mb of sequencing data\n";
	
	my $cortexFile = getCortexPathForSampleID($prefix, $sampleID, $k);
	my $cortexFile_flag = getCortexPathForOKFlag($prefix, $sampleID, $k);
	unless(-e $cortexFile_flag)
	{
		my $pe1 = $cortexFile . '.pe1';
		my $pe2 = $cortexFile . '.pe2';
		my $coverage = $cortexFile . '.coverage';
		my $coverageT = $cortexFile . '.coverageT';
	
		unlink($cortexFile) if (-e $cortexFile);
		unlink($coverage) if (-e $coverage);
		
		if(-e $coverageT)
		{
			unlink($coverageT) or die "Cannot unlink $coverageT";
		}
		
		open(PE1, '>', $pe1) or die "Cannot open $pe1";
		open(PE2, '>', $pe2) or die "Cannot open $pe2";
		foreach my $FASTQpair (@$FASTQs_aref)
		{
			die unless(scalar(@$FASTQpair) == 2);
			print PE1 $FASTQpair->[0], "\n";
			print PE2 $FASTQpair->[1], "\n";
		}
		
		close(PE1);
		close(PE2);
		my $cortex_bin = get_cortex_bin_by_nColours(1);
		my $cortex_cmd = qq($cortex_bin --mem_height $oneSample_cortex_height --mem_width 100 --pe_list $pe1,$pe2 --dump_binary $cortexFile --kmer_size $k --remove_pcr_duplicates --quality_score_threshold 5 --dump_covg_distribution $coverage);
		system($cortex_cmd) and die "Command $cortex_cmd failed";
		unlink($pe1);
		unlink($pe2);
		
		my $plot_cmd = qq(Rscript plotCoverage.R $coverage $sampleID);
		system($plot_cmd) and die "Command $plot_cmd failed";
		
		die "File $coverageT was not produced" unless(-e $coverageT);
		open(T, '<', $coverageT) or die "Cannot open $coverageT";
		my $T_line = <T>;
		chomp($T_line);
		close(T);
		die "File $coverageT has unexpected content, should have one integer" unless($T_line =~ /^\d+$/);
		
		system("echo 1 > $cortexFile_flag") and die "Echo command failed";
	}
	
	return $totalBases;
}

sub getCortexPathForOKFlag
{
	my $prefix = shift;
	my $sampleID = shift;
	my $k = shift;
	die unless(defined $k);
	my $cortexFile = getCortexPathForSampleID($prefix, $sampleID, $k);
	my $cortexFile_flag = $cortexFile . '.OK';
	return $cortexFile_flag;
}

sub getCortexPathForSampleID
{
	my $prefix = shift;
	my $sampleID = shift;
	my $k = shift;
	my $checkOK = shift;
	die unless(defined $k);
	
	if($checkOK)
	{
		my $flag_fn = getCortexPathForOKFlag($prefix, $sampleID, $k);
		die "Flag $flag_fn for sample $sampleID not there, abort." unless(-e $flag_fn);
	}
	
	return $cortex_temp_dir . "/${prefix}_${sampleID}_${k}.ctx";
}

sub getkMerThreshold
{
	my $prefix = shift;
	my $sampleID = shift;
	my $k = shift;	
	
	my $cortexFile = getCortexPathForSampleID($prefix, $sampleID, $k);
	my $coverageT = $cortexFile . '.coverageT';	
	
	open(T, '<', $coverageT) or die "Cannot open $coverageT";
	my $T_line = <T>;
	chomp($T_line);
	close(T);
	die "File $coverageT has unexpected content, should have one integer" unless($T_line =~ /^(\d+)$/);
	my $T = $1;
	
	return $T;
}

sub getLongestReadLength_FASTQ
{
	my $fn = shift;

	my $longestRead = 0;
	if($fn =~ /\.gz/)
	{
		die "Cannot open pipe to read output from zcat $fn";
	}
	
	open(F, '<', $fn) or die "Cannot open $fn";
	while(<F>)
	{
		chomp;
		next unless($_);
		die "Is file $fn really FASTQ? Weird beginning of line $.: ".substr($_, 0, 1) unless(substr($_, 0, 1) eq '@');
		my $S = <F>;
		chomp($S);
		my $L = length($S);
		my $p = <F>;
		die unless(substr($p, 0, 1) eq '+');
		<F>;
		
		$longestRead = $L if($L > $longestRead);
	}
	close(F);
	die unless($longestRead);
	
	return $longestRead;
}

sub countBasesInFASTQ
{
	my $fn = shift;
	my $bases = 0;
	if($fn =~ /\.gz/)
	{
		open(F, '-|', "zcat $fn") or die "Cannot open pipe to read output from zcat $fn";
	}
	else
	{
		open(F, '<', $fn) or die "Cannot open $fn";
	}
	
	while(<F>)
	{
		chomp;
		next unless($_);
		die "Is file $fn really FASTQ? Weird beginning of line $.: ".substr($_, 0, 1) unless(substr($_, 0, 1) eq '@');
		my $S = <F>;
		chomp($S);
		$bases += length($S);
		my $p = <F>;
		die unless(substr($p, 0, 1) eq '+');
		<F>;
	}
	close(F);
	return $bases;
}


sub countReadsInFASTQ
{
	my $fn = shift;
	my $reads = 0;
	if($fn =~ /\.gz/)
	{
		open(F, '-|', "zcat $fn") or die "Cannot open pipe to read output from zcat $fn";
	}
	else
	{
		open(F, '<', $fn) or die "Cannot open $fn";
	}
	
	while(<F>)
	{
		chomp;
		next unless($_);
		die "Is file $fn really FASTQ? Weird beginning of line $.: ".substr($_, 0, 1) unless(substr($_, 0, 1) eq '@');
		my $S = <F>;
		$reads++;
		my $p = <F>;
		die unless(substr($p, 0, 1) eq '+');
		<F>;
	}
	close(F);
	return $reads;
}

sub fastq_to_fasta
{
	my $input_fn = shift;
	my $output_fn = shift;
	die if($input_fn eq $output_fn);
	open(IN, '<', $input_fn) or die "Cannot open $input_fn";
	open(OUT, '>', $output_fn) or die "Cannot open $output_fn";
	while(<IN>)
	{
		chomp;
		next unless($_);
		die "Is file $input_fn really FASTQ? Weird beginning of line $.: ".substr($_, 0, 1) unless(substr($_, 0, 1) eq '@');
		my $readID = substr($_, 1);
		$readID =~ s/\s.+//g;
		
		my $S = <IN>;
		chomp($S);
		
		my $p = <IN>;
		die "Weird input in FASTQ file $input_fn line $.: expected a plus, got $p" unless(substr($p, 0, 1) eq '+');
		<IN>;
		
		print OUT '>', $readID, "\n", $S, "\n";
	}	
	close(IN);
	close(OUT);
}


sub observed_intersection_to_identity
{
	my $k = shift;
	my $total_kmers = shift;
	my $surviving_kMers = shift; 
	
	die unless($total_kmers > 0);
	die unless($surviving_kMers >= 0);
	die unless($surviving_kMers <= $total_kmers);
	
	my $observed_p = $surviving_kMers / $total_kmers;
	
	my $identity = $observed_p**(1/$k);
}

sub observed_kmer_intersection_likelihood
{
	my $k = shift;
	my $total_kMers = shift;
	my $identity = shift; 
	my $surviving_kMers = shift;
	
	die unless($total_kMers > 0);
	die unless(($identity >= 0) and ($identity <= 1));	
	die unless($surviving_kMers >= 0);
	die unless($surviving_kMers <= $total_kMers);
	
	my $kMer_survival_rate = $identity ** $k;

	# gsl_ran_binomial_pdf($k, $p, $n)
	my $p = gsl_ran_binomial_pdf($surviving_kMers, $kMer_survival_rate, $total_kMers);
	die unless($p <= 1);
	if($p == 0)
	{
		return -10000000;
	}
	else
	{
		return log($p);
	}
}	

sub get_cortex_bin_by_nColours
{
	my $ncolours = shift;
	die unless(defined $ncolours);
	if($ncolours <= 20)
	{
		return "${cortex_bin_dir}/cortex_var_31_c20";
	}
	elsif($ncolours <= 40)
	{
		return "${cortex_bin_dir}/cortex_var_31_c40";
	}
	elsif($ncolours <= 60)
	{
		return "${cortex_bin_dir}/cortex_var_31_c60";
	}	
	else
	{
		die "Requested $ncolours colors - don't have support for so many colours";
	}
}
