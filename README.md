# UltraPlexer


## Introduction
Sequencing by Oxford Nanopore technology is not only relative expensive, but also limited to a multiplex of 24 isolates, due t available barcodes.

To increasse the number of isolates that can be sequence simultaniously and thereby decreasing the per isolate sequencing costs we developed the UltraPlexer. A tool that matches non-barcoded long-read data of Oxford Nanopore Technologies to barcoded short-tead data of Illumina Technologies, based on k-mer frequencies, and assigns them to the matching isolates. The classifying algorithm has an error-rate of nearly 0% for a multi-species sequencing-run and an error-rate of below 1% for a multi-strain sequencing-run.


## Programm Requirements and Installation
**The following softwares need to be installed:**

- Perl
- R
- Cortex

At the moment the pathways to cortex still need to beintegrated manually into the UltraPlexer script, before running (UltraPlexer.pl, line 1412-1420). 


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

The programm itself, called by perl.

**--prefix prefix1**

Your chosen name for the UltraPlexer-run.

**--action classify / generateCallFile**

The command to classify the long-reads (classify) or to generate an output file from the classified data (generateCallFile).

**--samples_file /path/to/samplefile/samplefile1.txt**

A file containing the isolate ID, the pathway to the illumina_R1.fastq file and the pathway to the illumina_R2.fastq file. One line per isolate. The three infos seperated by tabs.

Example:
```
Isolate_1	/Data/isolate_1_R1.fastq	/Data/isolate_1_R2.fastq
Isolate_17	/Data/isolate_17_R1.fastq	/Data/isolate_17_R2.fastq
MRSA_H4	        /Data/MRSA_H4_R1.fastq          /Data/MRSA_H4_R2.fastq
Benjamin	/Data/Benjamin_R1.fastq	        /Data/Benjamin_R2.fastq
…
```

**--longReads_FASTQ /path/to/longreads/longreads1.fastq**

A file containing all long-reads from isolates in the “samplefile1.txt” that should be classified in standard fastq format.

**--classificationSource random**

The command to randomly distribute the long-reads, instead of classifying them.


## Output

**mixed_bacteria_10x.classification_k19.done**

This file is produced, when the UltraPlexer finished running correctly.

**mixed_bacteria_10x.classification_k19**

This file is produced by the first programmcall and contains all data produced while classifying the long-reads.

**mixed_bacteria_10x.classification_k19.called_kmers**

This file produces by the second programmcall contains the header information of the classified reads, the ID of the isolate the read is assigned to and the propability (?) for the by UltraPlexer classified reads. Onel ine per read. The three infos seperated by tabs.

Example: 
```
Read_1	        MRSA_H4         0.886902934417435
Read_332	Isolate_1	0.906668691485839
Read_336	Benjamin	0.895056000007794
Read_4100	Isolate_1	0.912532884787109
…
```

**mixed_bacteria_10x.classification_k19.called_random**

This file produces by the third programmcall contains the header information of the classified reads, the ID of the isolate the read is assigned to and the propability (?) for the random classified reads. Onel ine per read. The three infos seperated by tabs.

Example: 
```
Read_1          Benjamin	0.886902934417435
Read_332	Isolate_17	0.906668691485839
Read_336	Isolate 17 	0.895056000007794
Read_4100	MRSA_H4         0.912532884787109
…
```
