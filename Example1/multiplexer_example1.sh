#!/bin/bash
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=71:59:00
#PBS -A "NGS-Staphylo"
#PBS -N "multi_ex1"

set -e
 
cd /gpfs/project/dilthey/projects/multiplexer/
export PATH=$PATH:/gpfs/project/dilthey/software/racon/bin/:/gpfs/project/dilthey/software/wgsim:/gpfs/project/dilthey/software/PBSIM-PacBio-Simulator/src:/gpfs/project/dilthey/software/Unicycler:/gpfs/project/dilthey/software/MUMmer3.23

module load Perl
module load R

/usr/bin/time -v perl multiplexer.pl --prefix example1_plasmid --action classify --samples_file /gpfs/project/dilthey/projects/multiplexer_pool/Ultraplexing_example/example1_plasmid_ids_and_pathways.txt 	--longReads_FASTQ /gpfs/project/dilthey/projects/multiplexer_pool/Ultraplexing_example/example1_plasmid_read_pool.fastq
/usr/bin/time -v perl multiplexer.pl --prefix example1_plasmid --action generateCallFile --samples_file /gpfs/project/dilthey/projects/multiplexer_pool/Ultraplexing_example/example1_plasmid_ids_and_pathways.txt

