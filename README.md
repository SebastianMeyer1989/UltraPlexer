# UltraPlexer


## Introduction

UltraPlexing is a highly effective method for multiplexed long-read sequencing in the context of hybrid microbial genome assembly. Ultraplexing removes the need for molecular barcodes and assigns each long read to the short-read assembly graph it is most compatible with. While maintaining excellent assembly quality, Ultraplexing enables at least a doubling of the maximum number of samples per flow cell on the Nanopore platform, and a reduction in reagent costs and hands-on-time by a factor of 2.

To apply the UltraPlexing approach, simply pool equal amounts of DNA from the samples you want to multiplex, generate long-read sequencing data, and use the UltraPlexing algorithm to demultiplex the data. In order to apply UltraPlexing, short-read sequencing data for the same samples needs to be available at the time of analysis. If possible, optimize the DNA extration and library preparation processes for read length, as the ability of the UltraPlexing algorithm to assign long reads to isolates improves with increasing read lengths.

A preprint with accuracy evaluations will be made available soon.

## Program Requirements and Installation
The programm was tested on the following Operating System: CentOS Linux release 7.5.1804 (Core).

### The following programming languages and packages need to be installed:

- R
- [cortex_var](http://cortexassembler.sourceforge.net/index_cortex_var.html)
- Perl
- Perl Modules:
  - List::MoreUtils qw/mesh all uniq/
  - Math::GSL::Randist qw/gsl_ran_binomial_pdf/

Please modify `UltraPlexer.pl` so that it contains the correct path to your installation of Cortex (line 16). The algorithm expects to find Cortex binaries for k = 31 with 20, 40 and 60 colors.

## Running the UltraPlexer

### 1. Classify long-reads:
```
perl UltraPlexer.pl --prefix prefix1 --action classify --samples_file /path/to/samplefile/samplefile1.txt --longReads_FASTQ /path/to/longreads/longreads1.fastq
```
### 2. Generate human-readable classified output from classification file:
```
perl UltraPlexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt
```
### 3. Generate human-readable random output from classification file (as a comparison):
```
perl UltraPlexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt --classificationSource random
```

## Input

**perl UltraPlexer.pl**

The UltraPlexing algorithm.

**--prefix prefix1**

Your chosen prefix for the UltraPlexer run.

**--action classify / generateCallFile**

The command to classify the long-reads (`classify`) or to generate an output file from the classified data (`generateCallFile`).

**--samples_file /path/to/samplefile/samplefile1.txt**

A tab-separated file containing the isolate ID, the path to the illumina_R1.fastq file, and the path to the illumina_R2.fastq file. One line per isolate.

#### Example:
```
Isolate_1	/Data/isolate_1_R1.fastq	/Data/isolate_1_R2.fastq
Isolate_17	/Data/isolate_17_R1.fastq	/Data/isolate_17_R2.fastq
MRSA_H4	        /Data/MRSA_H4_R1.fastq          /Data/MRSA_H4_R2.fastq
Benjamin	/Data/Benjamin_R1.fastq	        /Data/Benjamin_R2.fastq
…
```

**--longReads_FASTQ /path/to/longreads/longreads1.fastq**

A FASTQ file containing the long reads to be classified.

**--classificationSource random**

The command to generate a random assignment of long reads to isolates (useful for benchmarking).

## Output

(examplary for prefix `mixed_bacteria_10x`):

**mixed_bacteria_10x.classification_k19.done**

This flag file is produced when the UltraPlexer finished running correctly.

**mixed_bacteria_10x.classification_k19**

This file is produced by the `classify` command and contains intermediate read classification data.

**mixed_bacteria_10x.classification_k19.called_kmers**

This file is produced by the `generateCallFile` command and contains, for each read, the isolate it has been assigned to, and a quality metric.

#### Example: 
```
Read_1	        MRSA_H4         0.886902934417435
Read_332	Isolate_1	0.906668691485839
Read_336	Benjamin	0.895056000007794
Read_4100	Isolate_1	0.912532884787109
…
```

**mixed_bacteria_10x.classification_k19.called_random**

This file is produced when specifying the `--classificationSource random` option. It contains a random allocation of reads to isolates.

## Creating fastq-files for further hybrid assemblies

(exemplary for prefix `mixed_bacteria_10x`)
 
### Call:
```
perl create_kmer_based_fastq_for_real_data.pl mixed_bacteria_10x.classification_k19.called_kmers path/to/longreads/longreads1.fastq  mixed_bacteria_10x
```

### Input:

**perl create_kmer_based_fastq_for_real_data.pl**

The script that produces fastq files from the calling table.

**mixed_bacteria_10x.classification_k19.called_kmers**

The calling table from the UltraPlexer run.

**path/to/longreads/longreads1.fastq**

The used long-read file.

**mixed_bacteria_10x**

Prefix for the run.

### Output:

**mixed_bacteria_10x-isolate1-predicted_reads.fastq**

A fastq file named after the run (mixed_bacteria_10x) and the isolate ID (isolate1), ending with “predicted_reads.fastq”.

## Example Run

In the following we exemplary describe how to simulate reads, create a read-pool, run the Ultraplexing algorithm and assemble the assigned reads, on basis of three random plasmids.
The nessecary data can be found in the folder “Example1”: 
- Three fasta files (Plasmid1.fna, Plasmid2.fna and Plasmid3.fna)
- A list of these three fasta files (example1_list_of_plasmids.txt)
- Different scripts (.pl and .R), needed to run the whole pipeline

After downloading the Folder "Example1" you just need to switch to it via terminal and call the necessary commands in the further  explained order (assumed, that all the requirements are met).


If you don't want to simulate data or the simulation is just not possible on your device, you can skip Step 1 and 2 (simulation and creation of the long-read pool) and use the long-read pool and other needed data we provided.
Therefor you need to download the zip-file "`Example1_simulated_example_data.zip`" from

https://uni-duesseldorf.sciebo.de/s/oHFl3FCArhPhHb5

and copy the content of the unzipped Folder "`Example1_simulated_example_data`" (Sim_Pipeline, example1_plasmid_ids_and_pathways.txt, example1_plasmid_read_pool.fastq and example1_plasmid_stats.txt) into your Folder "`Example1`" and continue the pipeline from step 3. "Ultraplexing:...".

### 1. Simulation (59s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl
- pbsim
- wgsim

#### Call:
```
perl simulation_pipeline.pl example1_list_of_plasmids.txt 8500 150
```
Important: Pathways for pbsim, pbsim qc-model and wgsim need to be replaced in the script "simulation_pipeline.pl" (line 10-12) to fit your installations, before running it.

#### Input:
- Genome files to simulate reads from
- List with the filenames of these genomes (example1_list_of_plasmids.txt)
- Desired mean read length of the long reads (8500)
- Desired coverage of the long reads (150)

Important: Please keep these numerical parameters for the example run as they are, since the scripts are still hard-coded at the moment for this mean lengh and coverage.

#### Output:
- A new folder (Sim_Pipeline) containing
  - Folder with data for each simulated genome (Plasmid1_l8500_c150 for example)
  - File with genomes files, the program could not find (missing_files.log)
  - File with the used simulation parameters (simulation_parameters.log)
- The final simulated reads are found in the files containing the tag “filtered”. For the genome “Plasmid1” they would be called:
  - “Plasmid1_l8500_c150-filtered_R1.fastq” (Illumina R1 short-reads)
  - “Plasmid1_l8500_c150-filtered_R1.fastq” (Illumina R2 short-reads)
  - “Plasmid1_l8500_c150-filtered_complete.fastq” (Nanopore long-reads)

### 2. Creating a shuffled long-read pool (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl

#### Call:
```
perl create_pool.pl example1_list_of_plasmids.txt 3 10000000 example1_plasmid
```
#### Input:
- List with the filenames of these genomes (example1_list_of_plasmids.txt)
- Number of genomes that should be pooled together from that list (3)
- Number of bases that should be pooled together in total (10000000)
- Prefix for the run (example1_plasmid)

#### Output:
- Shuffled long-read pool (example1_plasmid_read_pool.fastq)
- File containing parameters of the pooling step (example1_plasmid_stats.txt)
- File containing the IDs and simulated short-read pathways of all pooled genomes, needed for the next step (example1_plasmid_ids_and_pathways.txt)

### 3. Ultraplexing: Classification of long-reads (2m32s runtime, 1 CPUs, <7gb used memory)

#### Requirements:
- perl
- R

#### Call:
```
perl UltraPlexer.pl --prefix example1_plasmid --action classify --samples_file example1_plasmid_ids_and_pathways.txt --longReads_FASTQ example1_plasmid_read_pool.fastq
```
#### Input:
- Prefix for the run (example1_plasmid)
- Action the script should do (classify)
- Tab seperated file containing the IDs and absolute pathways of the simulated short-read of all pooled genomes, one genome per line (example1_plasmid_ids_and_pathways.txt). This was produced by the pooling script.
- Shuffled long-read pool (example1_plasmid_read_pool.fastq). This was produced by the pooling script.

#### Output:
- File with classified long-reads (example1_plasmid.classification_k19)
- File produced when the classification step finished correctly (example1_plasmid.classification_k19.done)
- Folder containing temporary files produced by cortex (cortex_temp)

### 4. Ultraplexing: Creating the assignment table (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl
- R

#### Call:
```
perl UltraPlexer.pl --prefix example1_plasmid --action generateCallFile --samples_file example1_plasmid_ids_and_pathways.txt
```
#### Input:
- Prefix for the run (example1_plasmid)
- Action the script should do (generateCallFile)
- Tab seperated file containing the IDs and absolute pathways of the simulated short-read of all pooled genomes, one genome per line (example1_plasmid_ids_and_pathways.txt). This was produced by the pooling script.

#### Output:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)


### 5. Creating fastq-files for each genome (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl

#### Call:
```
perl create_kmer_based_fastq_for_simulations.pl example1_plasmid.classification_k19.called_kmers example1_plasmid_read_pool.fastq example1_plasmid
```
#### Input:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)
- Shuffled long-read pool (example1_plasmid_read_pool.fastq)
- Prefix for the run (example1_plasmid)

#### Output:
- Two fastq-files for each simulated genome:
  - Predicted reads for the genome in a file, ending with “...predicted_reads.fastq”.
  - True reads for the genome in a file, ending with “...true_reads.fastq”. These are the reads, that were originally simulated for the genome.

### 6. Comparing true and predicted reads (1s runtime, 1 CPUs, <1gb used memory)

(examplary for `Plasmid1`)

#### Requirements:
- perl

#### Call:
```
perl parse_calling_tbl.pl  example1_plasmid.classification_k19.called_kmers
```
#### Input:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)

#### Output:
- File containing primary stats (example1_plasmid.classification_k19.called_kmers.stats). These are not necessarily important for you.
- File containing final stats (example1_plasmid.classification_k19.called_kmers.stats2). This file shows the name of the used calling table and the number and ratio of correct called reads in total (first line) and the number of falsly and correctly called reads for every used genome (following lines).
##### Example:
```
example1_plasmid.classification_k19.called_kmers Summary		Correct_Reads:	1147	False_Reads:	47	Ratio_Correct_Reads:	0.960636515912898
Plasmid1	false: 4		Plasmid1	true: 394
Plasmid2	false: 16		Plasmid2	true: 382
Plasmid3	false: 27		Plasmid3	true: 371
```
Here you can see, that over 96% of the simulated reads were assigned correctly. The missing <4% are miss-assignments due to sequence homology, or do not really affect hybrid assemblies (as we found out in our experiments).

### 7. Hybrid assembly (1-2h runtime per assembly, 2 CPUs, <2gb used memory)

(examplary for predicted reads for `Plasmid1`)

#### Requirements:
- perl
- Python
- Unicycler
- SPAdes
- Racon
- Pilon
- Bowtie2
- SamTools
- BLASTp

#### Call:
```
/gpfs/project/dilthey/software/Unicycler/unicycler-runner.py --spades_path /software/SPAdes/3.11.1/ivybridge/bin/spades.py --racon_path /gpfs/project/dilthey/software/racon/bin/racon --pilon_path /software/pilon/1.22/pilon-1.22.jar -t 2 -1 Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R1.fastq -2 Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R2.fastq -l example1_plasmid-Plasmid1-predicted_reads.fastq -o example1_plasmid-Plasmid1-predicted_reads_unicycler
```
Important: Pathways for unicycler, spades, racon and pilon need to be replaced to fit your installations.

#### Input:
- Path to Unicycler (/gpfs/project/dilthey/software/Unicycler/unicycler-runner.py)
- Path to SPAdes (/software/SPAdes/3.11.1/ivybridge/bin/spades.py)
- Path to Racon (/gpfs/project/dilthey/software/racon/bin/racon)
- Path to Pilon (/software/pilon/1.22/pilon-1.22.jar)
- Number of threads used (2)
- Path to the Illumina R1 read file (Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R1.fastq)
- Path to the Illumina R2 read file (Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R2.fastq)
- Path to the long-read file (example1_plasmid-Plasmid1-predicted_reads.fastq)
- Prefix for the output folder (example1_plasmid-Plasmid1-predicted_reads_unicycler)

#### Output:
- Assembly-folder for each assembled genome, for example called “example1_plasmid-Plasmid1-predicted_reads_unicycler”, containing:
  - De-bruijn (Bandage) graphs for all assembly steps (“001_best_spades_graph.gfa” to “007_rotated.gfa” and “assembly.gfa”).
  - Assembly file in fasta format (assembly.fasta)
  - Log file with parameters of the assembly run (unicycler.log)

### 8. Comparing assemblies of true and predicted reads using nucmer (1s runtime, 1 CPUs, <1gb used memory)

(examplary for `Plasmid1`)

#### Requirements:
- nucmer
- delta-filter

#### Calls:
1.
```
nucmer -p Plasmid1 example1_plasmid-Plasmid1-predicted_reads_unicycler/assembly.fasta example1_plasmid-Plasmid1-true_reads_unicycler/assembly.fasta
```
2.
```
mummerplot2 Plasmid1.delta --png -p Plasmid1
```
#### Input:
1.
- Prefix for the run (Plasmid1)
- Predicted read assembly (example1_plasmid-Plasmid1-predicted_reads_unicycler/assembly.fasta)
- True read assembly  (example1_plasmid-Plasmid1-true_reads_unicycler/assembly.fasta)
2.
- Delta file from previous nucmer run (Plasmid1.delta)
- Prefix for the run (Plasmid1)

#### Output:
1.
- Alignment of true read assembly and predicted read assembly (Plasmid1.delta)
2.
- Graphical representation of the delta file (Plasmid1.png)

If you take a look at the produced graphics, you will see, that the assemblies of the predicted reads align perfectly to the assemblies of the true reads.
