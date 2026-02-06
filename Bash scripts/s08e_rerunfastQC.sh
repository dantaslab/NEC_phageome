#!/bin/bash
#===============================================================================
#
# File Name    : s04_rerunfastqc.sh
# Description  : This script will run fastqc in parallel
# Usage        : sbatch s04_rerunfastqc.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2020-04-28 by Bejan Mahmud
# Last Modified: 2022-06-14
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=rerunfastqc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-2120
#SBATCH --output=slurm/fastQC/x_rerunfastqc_%a.out
#SBATCH --error=slurm/fastQC/y_rerunfastqc_%a.err

module load fastqc

#basedir
basedir="/scratch/gdlab/kailun/Phageome/NEC/Reads_mapping"
indir="${basedir}/NEC_cleanreads_time"
outdir="${basedir}/fastqc"

#make directory
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/sample_list.txt`

#running fastqc
fastqc -o ${outdir} -f fastq ${indir}/${sample}_R1_CLEAN.fastq.gz
fastqc -o ${outdir} -f fastq ${indir}/${sample}_R2_CLEAN.fastq.gz
