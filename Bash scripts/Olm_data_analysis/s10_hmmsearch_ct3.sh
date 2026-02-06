#!/bin/bash
#===============================================================================
# File Name    : s10_hmmsearch_ct3.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s10_hmmsearch_ct3.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Jul 22 2021 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2026-09-13
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=HMMER
#SBATCH --mem=640G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm/hmmsearch/x_hmmsearch_%A_%a.out
#SBATCH --error=slurm/hmmsearch/y_hmmsearch_%A_%a.out

#module load hmmer
eval $( spack load --sh hmmer )

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d10_hmmsearch_ct3"


hmmsearch --tblout ${basedir}/d10_hmmsearch_ct3/hmmsearch_SRR_ct3.out -E 0.05 ${indir}/final_list.hmms ${basedir}/d08_CheckV_ct3/tmp/proteins.faa