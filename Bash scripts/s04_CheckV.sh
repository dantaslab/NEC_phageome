#!/bin/bash
#===============================================================================
# File Name    : s04_CheckV.sh
# Description  : Assessing the quality of metagenome-assembled viral genomes
# Usage        : sbatch s04_CheckV.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Modified     : Nov 02 2021
# Created      : May 14 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkv
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/checkv/x_checkv_%A.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/checkv/y_checkv_%A.out
#SBATCH -C cpu_E52650

module load miniconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /scratch/ref/gdlab/CheckV
export CHECKVDB=/scratch/ref/gdlab/CheckV/checkv-db-v1.0


set -x
checkv end_to_end /scratch/gdlab/kailun/Phageome/NEC/VAMB/Concate_bins/Pass4_viral_concate_contig2.fna /scratch/gdlab/kailun/Phageome/NEC/CheckV_concate_pass4/ -t 16

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
