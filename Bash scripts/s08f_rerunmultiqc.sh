#!/bin/bash
#===============================================================================
#
# File Name    : s02_rerunmultiqc.sh
# Description  : This script will run fastqc in parallel
# Usage        : sbatch s05_multiqc.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2020-04-28 by Bejan Mahmud
# Last Modified: 2022-06-14
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=rmultiqc
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --output=slurm/multiqc/x_multiqc_%a.out
#SBATCH --error=slurm/multiqc/y_multiqc_%a.err

module load multiqc

#basedir
basedir="/scratch/gdlab/kailun/Phageome/NEC/Reads_mapping/fastqc"

#running multiqc
multiqc ${basedir}/