#!/bin/bash
#===============================================================================
# File Name    : s11a_busco.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s11a_busco.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Oct 04 2021
# Last Modified: 2024-09-28
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=busco
#SBATCH --array=1-8674%500
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --output=slurm/busco/z_busco_%A_%a.out
#SBATCH --error=slurm/busco/z_busco_%A_%a.out

#activate env
#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#conda activate busco_v5.2.2

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm/BUSCO"
indir="${basedir}/pass3_single_contigs"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240928_Olm_srr_phage_pass3_list_ct3.txt`


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
