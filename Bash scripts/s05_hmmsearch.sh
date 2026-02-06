#!/bin/bash
#===============================================================================
# File Name    : s05_hmmsearch.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s05_hmmsearch.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Jul 22 2021 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: Nov 03 2021
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=HMMER
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/hmmsearch/x_hmmsearch_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/hmmsearch/y_hmmsearch_%A_%a.out

module load hmmer

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/HmmSearch"


hmmsearch --tblout ${indir}/hmmsearch.out -E 0.05 ${indir}/final_list.hmms ${basedir}/CheckV_concate_pass4/tmp/proteins.faa
