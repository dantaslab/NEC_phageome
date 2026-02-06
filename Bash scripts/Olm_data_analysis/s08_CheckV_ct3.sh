#!/bin/bash
#===============================================================================
# File Name    : s08_CheckV_ct3.sh
# Description  : Assessing the quality of metagenome-assembled viral genomes
# Usage        : sbatch s08_CheckV_ct3.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 3.0
# Modified     : Nov 02 2021
# Created      : 2024-09-13
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkv
#SBATCH --cpus-per-task=16
#SBATCH --mem=320G
#SBATCH --output=slurm/checkv/x_checkv_ct3_SRS_%A.out
#SBATCH --error=slurm/checkv/y_checkv_ct3_SRS_%A.out

eval $( spack load --sh /buvzc6u )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /ref/gdlab/software/envs/CheckV
export CHECKVDB=/ref/gdlab/software/envs/CheckV/checkv-db-v1.5


set -x
checkv end_to_end /scratch/gdlab/kailun/Phageome/NEC_Olm/d07_cenotetaker3/CT3_fasta/240913_Olm_SRR_CT3.fasta /scratch/gdlab/kailun/Phageome/NEC_Olm/d08_CheckV_ct3/ -t 16

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
