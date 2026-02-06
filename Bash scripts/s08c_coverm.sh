#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm_v3.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12c_coverm_v3.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 3.0
# Created On   : 2022-06-01
# Last Modified: 2022-06-14
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=coverm
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-144
#SBATCH --output=slurm/coverm/x_coverm_%A_%a.out
#SBATCH --error=slurm/coverm/y_coverm_%A_%a.out

module load coverm

basedir="/scratch/gdlab/kailun/Phageome/NEC"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_144_new.txt`
#ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie_unbin/ref_list.txt`
indir="${basedir}/Reads_mapping"
outdir="${indir}/out_files_tpmean"

#make output directory
mkdir -p ${outdir}
# relative_abundance mean trimmed_mean covered_bases covered_fraction variance length count reads_per_base
time coverm genome --bam-files ${indir}/mapping_files/${sample}_bam_files/*.bam --genome-fasta-files ${indir}/phage_contigs/${sample}_contigs/*.fna -m relative_abundance mean trimmed_mean covered_bases covered_fraction variance length count reads_per_base --min-read-percent-identity 95 --min-covered-fraction 0 -o ${outdir}/${sample}_coverm_output.tsv 

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
