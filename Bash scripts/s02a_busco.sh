#!/bin/bash
#===============================================================================
# File Name    : s01_busco.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s01_busco.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Oct 04 2021
# Last Modified: May 24 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=busco
#SBATCH --array=1-11094
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/busco/x_busco_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/busco/y_busco_%A_%a.out

#activate env
#conda activate /scratch/ref/gdlab/BUSCO/5.2.2

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/VAMB/Concate_bins"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${indir}/concate_list.txt`


set -x
time busco -i ${indir}/${sample}.fasta -l bacteria_odb10 -o ${sample} -m genome -e 0.05
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi
