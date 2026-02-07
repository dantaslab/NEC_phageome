#!/bin/bash
#===============================================================================
#
# File Name    : s11_blast2seq.sh
# Description  : This script will run blastn
# Usage        : sbatch s11_blast2seq.sh
# Author       : Kailun Zhang
# Version      : 1.1
# Created On   : 2022-01-30
# Last Modified: 2024-04-21
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=blast
#SBATCH --array=1-144
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm/blast/x_blast2seq_%A_%a.out
#SBATCH --error=slurm/blast/y_blast2seq_%A_%a.out

eval $( spack load --sh /u7ssbm4 )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/mec/mec_genes"
outdir="${indir}"

#read in the slurm array task
sample1=`sed -n ${SLURM_ARRAY_TASK_ID}p ${indir}/mec_list1.txt`
sample2=`sed -n ${SLURM_ARRAY_TASK_ID}p ${indir}/mec_list2.txt`


set -x

time blastn -subject ${indir}/${sample1}.fasta -query ${indir}/${sample2}.fasta -perc_identity 0 -evalue 1000 -outfmt 6 -out ${outdir}/${sample1}_${sample2}_shortblastn.out

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi