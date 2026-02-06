#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12c_coverm.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Last Modified: 2022-06-01
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=coverm
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2120
#SBATCH --output=slurm/coverm/x_coverm_%A_%a.out
#SBATCH --error=slurm/coverm/y_coverm_%A_%a.out

module load coverm

basedir="/scratch/gdlab/kailun/Phageome/NEC"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/sample_list.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/ref_list.txt`
indir="${basedir}/d12_bowtie/d12_bowtie2/${sample}"
outdir="${basedir}/d12_bowtie/d12_coverm"

#make output directory
mkdir -p ${outdir}

time coverm genome --bam-files ${indir}/${sample}_sorted.bam --genome-fasta-files ${basedir}/d12_bowtie/${ref}_contigs.fasta -m relative_abundance mean trimmed_mean covered_bases --min-read-percent-identity 95 -o ${outdir}/${sample}_coverm_output.tsv 

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
