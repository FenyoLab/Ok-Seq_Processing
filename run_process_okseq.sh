#!/bin/bash
#SBATCH --mem=50gb                     # Job memory request
#SBATCH --time=48:00:00                # Time limit d-h:m:s  # should be 48 hrs
#SBATCH -p cpu_medium

# Example shell script for running python script on HPC cluster
process_dir=$1 # dir where fastq files are located; output of intermediate processing steps will be here
final_output_dir=$2 # dir to save the final output file and the log 
fastq_r1=$3 # fastq r1 file
fastq_r2=$4 # fastq r2 file, *set this to "." if it is SE*
samfile=$5 # sam file (should be .sam)
logfile=$6 # log file (should be .txt)
outfile=$7 # final output file (should be .txt)

module purge
module unload python
module load python/cpu/2.7.15-ES
module load samtools/1.10
module load picard/2.18.11
module load bedtools/2.27.1
module load deeptools/3.2.1
module load fastqc
module load trimgalore/0.5.0
module load bowtie2/2.4.1

python process_okseq_SE_PE.py $process_dir $final_output_dir $fastq_r1 $fastq_r2 $samfile $logfile $outfile 0

