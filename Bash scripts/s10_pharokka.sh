#!/bin/bash
#===============================================================================
# File Name    : s09_pharokka.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s09_pharokka.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 1.1 (v1.3.0)
# Modified     : 2023-07-09
# Created      : 2022-04-28
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=prokka
#SBATCH --array=1-3496%50
#SBATCH --mem=4G
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm/pharokka/x_pharokka_%A_%a.out
#SBATCH --error=slurm/pharokka/y_pharokka_%A_%a.err

#eval $( spack load --sh miniconda3 )
#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#conda activate pharokka_env_1.3.0

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/Final_viral_single_contig"
outdir="${basedir}/Pharokka_final_3469"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/230524_all_phage_3469_list.txt`

set -x
time python /ref/gdlab/software/envs/pharokka/v1.3.0/pharokka/bin/pharokka.py -i ${indir}/${sample}.fasta -o ${outdir}/${sample} -d /ref/gdlab/software/envs/pharokka/v1.3.0/pharokka/databases -t ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi