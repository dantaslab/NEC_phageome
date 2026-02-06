#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12c_coverm.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Last Modified: 2024-11-18
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=coverm
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=799,875
#SBATCH --output=slurm/coverm/z_coverm_%A_%a.out
#SBATCH --error=slurm/coverm/z_coverm_%A_%a.out

source /ref/gdlab/software/minicondaPhables/bin/activate
conda activate /ref/gdlab/software/minicondaPhables/envs/coverm

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/sample_list.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/ref_list.txt`
indir="${basedir}/d12_bowtie/d12_bowtie2/${sample}"
outdir="${basedir}/d12_bowtie/CoverM_contig"

#make output directory
mkdir -p ${outdir}

time coverm contig --bam-files ${indir}/${sample}_sorted.bam -m mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm --min-read-percent-identity 95 -o ${outdir}/${sample}_coverm_output.tsv 

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
