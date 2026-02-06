#!/bin/bash
#===============================================================================
# File Name    : s05_vibrant.sh
# Description  : Virus Identification By iteRative ANnoTation
# Usage        : sbatch s05_vibrant.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 20 2021
# Last Modified: May 14 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vibrant
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vibrant/x_vibrant_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vibrant/y_vibrant_%A_%a.out


module use /opt/htcf/modules
module use /opt/htcf/modules-legacy
module use /opt/apps/labs/gdlab/modules

module load miniconda3
CONDA_BASE=$(conda info --base)
source /home/kailun/conda/miniconda3/etc/profile.d/conda.sh

conda activate /scratch/ref/gdlab/VIBRANT/VIBRANT-master
module load hmmer/3.3.2
module load prodigal/2.6.3

python3.7 /scratch/ref/gdlab/VIBRANT/VIBRANT-master/VIBRANT/VIBRANT_run.py -i /scratch/gdlab/kailun/Phageome/NEC/VAMB/Concate_bins/all_concate_contigs.fna -t 20 -virome -folder /scratch/gdlab/kailun/Phageome/NEC/VIBRANT_concate
