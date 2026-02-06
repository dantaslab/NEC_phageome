#!/bin/bash
#===============================================================================
# File Name    : s02c_Concatenate.sh
# Description  : Concatenate the contigs
# Usage        : sbatch s02c_Concatenate.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.1
# Modified     : 2022-05-14
# Created      : 2022-05-24
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=Concatenate
#SBATCH --cpus-per-task=1
#SBATCH --array=1-10711
#SBATCH --mem=2G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/x_Concatenate_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/y_Concatenate_%A_%a.out

basedir="/scratch/gdlab/kailun/Phageome/NEC"
#indir="${basedir}/VAMB/All_bins"
indir="${basedir}/VAMB/vamb_out_end_1/bins"
outdir="${basedir}/VAMB/Concate_bins"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${indir}/seq_list.txt`

set +x

seqname=`grep -m1 "" ${indir}/${sample}.fna`
echo "${seqname}" > ${outdir}/${seqname}.fasta
grep -v '>' ${indir}/${sample}.fna >> ${outdir}/${seqname}.fasta


RC=$?
set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
