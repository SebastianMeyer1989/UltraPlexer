# UltraPlexer
## Requirements

- Perl
- R
- ?


## Programmcall

**Classification:**
```
perl multiplexer.pl --prefix prefix1 --action classify --samples_file /path/to/samplefile/samplefile1.txt --longReads_FASTQ /path/to/longreads/longreads1.fastq
```
**Generate classified output:**
```
perl multiplexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt
```
**Generate random output:**
```
perl multiplexer.pl --prefix prefix1 --action generateCallFile --samples_file /path/to/samplefile/samplefile1.txt --classificationSource random
```

## Explanation

Programcall 1 (classification) starts the classifying itself.

Programcall 2 (generate classified output) creates a human readable output file from the classified data, containing a list of the by UltraPlexer classified reads.

Programcall 3 (generate random output) creates a human readable output file from the classified data, containing a list of random classified reads (by the option “--classificationSource random”).


## Input

**multiplexer.pl**

The programm itself

**samplefile1.txt**

A file containing the isolate ID, the pathway to the illumina_R1.fastq file and the pathway to the illumina_R2.fastq file. One line per isolate. The three infos seperated by tabs.

Example:
```
Isolate_1	/Data/isolate_1_R1.fastq	/Data/isolate_1_R2.fastq
Isolate_17	/Data/isolate_17_R1.fastq	/Data/isolate_17_R2.fastq
MRSA_H4	        /Data/MRSA_H4_R1.fastq          /Data/MRSA_H4_R2.fastq
Benjamin	/Data/Benjamin_R1.fastq	        /Data/Benjamin_R2.fastq
…
```

**longreads1.fastq**

A file containing all long-reads from isolates in the “samplefile1.txt” that should be classified in standard fastq format.


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
