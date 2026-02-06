#!/bin/bash
#===============================================================================
# File Name    : s01.3_copyfile.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s01.3_copyfile.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 1.3
# Modified     : Aug 30 2021 
# Created      : 2024-09-13
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=copyfile
#SBATCH --mem=1G
#SBATCH --array=1-160
#SBATCH --output=slurm/copyfile/x_copyfile_%a.out
#SBATCH --error=slurm/copyfile/y_copyfile_%a.err

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d07_cenotetaker3"
outdir="${indir}/CT3_fasta"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

set -x
time cp ${indir}/${sample}/${sample}_virus_sequences.fna ${outdir}/${sample}_virus_sequences.fna
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi