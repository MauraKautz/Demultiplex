#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --time=1-0
#SBATCH --output=Part1R1_%j.out
#SBATCH --error=Part1R1_%j.err

/usr/bin/time -v ./demultiplex.py -f1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -f2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -f3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -f4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -q 20 -m ./matched_indexes.txt