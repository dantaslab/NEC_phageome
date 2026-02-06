#!/bin/bash
#===============================================================================
# File Name    : s05_vibrant.sh
# Description  : Virus Identification By iteRative ANnoTation
# Usage        : sbatch s05_vibrant.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 2.0
# Created On   : Nov 20 2021
# Last Modified: Jun 09 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vibrant
#SBATCH --cpus-per-task=20
#SBATCH --mem=70G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vibrant/x_vibrant_%A.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vibrant/y_vibrant_%A.out

# In a Linux Screen
#srun -c 20 --time=60:00:00 --mem=240G --pty bash -i 
eval $( spack load --sh python@3.9.12 )
eval $( spack load --sh hmmer )
eval $( spack load --sh prodigal )
eval $( spack load --sh /knknaya ) # loading py-seaborn 
eval $( spack load --sh /yl6lkpn ) # loading py-matplotlib

python3 /ref/gdlab/software/envs/VIBRANT-Master/VIBRANT/VIBRANT_run.py -i /scratch/gdlab/kailun/Phageome/NEC_Olm/d07_cenotetaker3/CT3_fasta/240913_Olm_SRR_CT3.fasta -t 20 -virome -folder /scratch/gdlab/kailun/Phageome/NEC_Olm/d09_VIBRANT_ct3

