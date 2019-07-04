# UltraPlexer


## Introduction

UltraPlexing is a highly effective method for multiplexed long-read sequencing in the context of hybrid microbial genome assembly. Ultraplexing removes the need for molecular barcodes and assigns each long read to the short-read assembly graph it is most compatible with. While maintaining excellent assembly quality, Ultraplexing enables at least a doubling of the maximum number of samples per flow cell on the Nanopore platform, and a reduction in reagent costs and hands-on-time by a factor of 2.

To apply the UltraPlexing approach, simply pool equal amounts of DNA from the samples you want to multiplex, generate long-read sequencing data, and use the UltraPlexing algorithm to demultiplex the data. In order to apply UltraPlexing, short-read sequencing data for the same samples needs to be available at the time of analysis. If possible, optimize the DNA extration and library preparation processes for read length, as the ability of the UltraPlexing algorithm to assign long reads to isolates improves with increasing read lengths.

A preprint with accuracy evaluations will be made available soon.

## Program Requirements and Installation
**The following programming languages and packages need to be installed:**

- Perl
- R
- [cortex_var](http://cortexassembler.sourceforge.net/index_cortex_var.html)

Please modify `UltraPlexer.pl` so that it contains the correct path to your installation of Cortex (line 16). The algorithm expects to find Cortex binaries for k = 31 with 20, 40 and 60 colors.

## Running the UltraPlexer

**1. Classify long-reads:**
```
perl multiplexer.pl --prefix prefix1 --action classify --samples_file /path/to/samplefile/samplefile1.txt --longReads_FASTQ /path/to/longreads/longreads1.fastq
```
**2. Generate human-readable classified output from classification file:**
```
perl multiplexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt
```
**3. Generate human-readable random output from classification file (as a comparison):**
```
perl multiplexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt --classificationSource random
```

## Input

**perl multiplexer.pl**

The UltraPlexing algorithm.

**--prefix prefix1**

Your chosen prefix for the UltraPlexer run.

**--action classify / generateCallFile**

The command to classify the long-reads (`classify`) or to generate an output file from the classified data (`generateCallFile`).

**--samples_file /path/to/samplefile/samplefile1.txt**

A tab-separated file containing the isolate ID, the path to the illumina_R1.fastq file, and the path to the illumina_R2.fastq file. One line per isolate.

Example:
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

(for prefix `mixed_bacteria_10x`):

**mixed_bacteria_10x.classification_k19.done**

This flag file is produced when the UltraPlexer finished running correctly.

**mixed_bacteria_10x.classification_k19**

This file is produced by the `classify` command and contains intermediate read classification data.

**mixed_bacteria_10x.classification_k19.called_kmers**

This file is produced by the `generateCallFile` command and contains, for each read, the isolate it has been assigned to, and a quality metric.

Example: 
```
Read_1	        MRSA_H4         0.886902934417435
Read_332	Isolate_1	0.906668691485839
Read_336	Benjamin	0.895056000007794
Read_4100	Isolate_1	0.912532884787109
…
```

**mixed_bacteria_10x.classification_k19.called_random**

This file is produced when specifying the `--classificationSource random` option. It contains a random allocation of reads to isolates.
