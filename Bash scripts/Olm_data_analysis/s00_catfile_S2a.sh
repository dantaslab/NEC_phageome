#!/bin/bash
#===============================================================================
# File Name    : s00_catfile.sh
# Description  : This script will concatenate read files
# Usage        : sbatch s00_catfile.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 1.0
# Modified     : 2024-06-18
# Created      : 2024-08-21
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=catfile
#SBATCH --mem=1G
#SBATCH --array=2-16
#SBATCH --output=slurm/catfile/a_catfile_S2_%a.out
#SBATCH --error=slurm/catfile/a_catfile_S2_%a.err

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/Rename/PRJNA376566"
outdir="${basedir}/SRA_Reads"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/SRA_Reads/240821_S2_infant_list.txt`

set -x
time 
cat ${indir}/${sample}*_1.fastq.gz > ${outdir}/${sample}_000G1_1.fastq.gz
cat ${indir}/${sample}*_2.fastq.gz > ${outdir}/${sample}_000G1_2.fastq.gz

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi