#!/bin/sh
#===============================================================================
#
# File Name    : s16c_humann3_join.sh
# Description  : Join the output files (gene families and abundance) from the HUMAnN runs from all samples
# Usage        : sbatch s16c_humann3_join.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Modified     : Kailun Zhang on 241202
# Version      : 1.1
# Created On   : Mon Jan  02 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=humann3
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/humann/m_humann3-merge_%a.out
#SBATCH --error=slurm/humann/m_humann3-merge_%a.out

#eval $( spack load --sh miniconda3@4.10.3 )
#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh 
#source activate /scratch/gdlab/gmark/conda_envs/humann3
source activate /ref/gdlab/software/envs/humann3

basedir="$PWD"
indir="${basedir}/d16_humann3/relab"


# Run Humann3 humann_renorm_table utility script
humann_join_tables --input ${indir} \
--output ${indir}/merged_humamm3_genefamilies_Olm.tsv --file_name genefamilies_relab


humann_join_tables --input ${indir} \
--output ${indir}/merged_humamm3_pathabundance_Olm.tsv --file_name pathabundance_relab