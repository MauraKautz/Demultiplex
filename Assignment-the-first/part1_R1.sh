#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --time=1-0
#SBATCH --output=Part1R1_%j.out
#SBATCH --error=Part1R1_%j.err

/usr/bin/time -v ./mean_qual.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l 101 -o R1_hist.png