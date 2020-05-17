#!/bin/bash
#SBATCH --mem=50gb                     # Job memory request
#SBATCH --time=48:00:00              # Time limit d-h:m:s  # should be 48 hrs, use cpu-medium when submitting batch job!

# Example shell script for running python script on HPC cluster

fastq_r1=$1 # full path to fastq r1 file
fastq_r2=$2 # full path to fastq r2 file, *set this to "." if it is SE*
samfile=$3 # full path to location of sam file
logfile=$4 # full path, or not if you want the log in current dir
outfile=$5 # full path, or not if you want the final output in current dir

module purge
module unload python
module load python/cpu/2.7.15-ES
module load samtools
module load picard
module load bedtools
module load deeptools
module load fastqc
module load trimgalore/0.5.0
module load bowtie2

python process_okseq_SE_PE.py $fastq_r1 $fastq_r2 $samfile $logfile $outfile 0

