#!/bin/bash
#===============================================================================
# File Name    : s07_vcontact2.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s07_vcontact2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Dec 14 2021
# Last Modified: June 02 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vcontact2
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vcontact2/x_vcontact2_%A.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vcontact2/y_vcontact2_%A.out

#module use /opt/htcf/modules
#module use /opt/htcf/modules-legacy
#module use /opt/apps/labs/gdlab/modules

module load miniconda3
CONDA_BASE=$(conda info --base)
source /home/kailun/conda/miniconda3/etc/profile.d/conda.sh

# activate env
conda activate /home/kailun/conda/vContact2

module load openjdk

basedir="/scratch/gdlab/kailun/Phageome/NEC"

set -x
vcontact2 --raw-proteins ${basedir}/concate_pass4_proteins.faa --rel-mode 'Diamond' --proteins-fp ${basedir}/concate_pass4_g2g.csv --db 'ProkaryoticViralRefSeq207-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/kailun/conda/MAVERICLab-vcontact2-34ae9c466982/bin/cluster_one-1.0.jar --output-dir ${basedir}/ConcatePass4_vConTACT_ref207_Outputs


RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in vContact2"
  exit $RC
fi