#!/bin/bash
#===============================================================================
# File Name    : s01_CenoteTaker2.sh
# Description  : Discover and annoate virus sequences/genomes
# Usage        : sbatch s01_CenoteTaker2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 2.0
# Modified     : Aug 17 2021 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Created      : May 05 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=ct2
#SBATCH --cpus-per-task=4
#SBATCH --array=1-5
#SBATCH --mem=30G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/CenoteTaker2/x_ct2_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/CenoteTaker2/y_ct2_%A_%a.err
#SBATCH -C cpu_E52650

module use /opt/htcf/modules
module use /opt/htcf/modules-legacy
module use /opt/apps/labs/gdlab/modules
module load miniconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /opt/apps/labs/gdlab/envs/cenote-taker2/2.1.3/cenote-taker2_env
CENOTE_BASE="/opt/apps/labs/gdlab/envs/cenote-taker2/2.1.3/Cenote-Taker2"

DATA_BASE="/scratch/gdlab/kailun/Phageome/NEC"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${DATA_BASE}/scaffold_list2.txt`
sample_out=`sed -n ${SLURM_ARRAY_TASK_ID}p ${DATA_BASE}/scaffold_list2_new.txt`
#mkdir /tmp/kailun_cenote

set -x
time ${CENOTE_BASE}/run_cenote-taker2.py \
    -c ${DATA_BASE}/metaspades/scaffolds_new/${sample}.fasta \
    -r ${sample_out} \
    -p True \
    -m 30 \
    -t 4
#    --scratch_directory /tmp/kailun_cenote
#    --known_strains blast_knowns \
#    --blastn_db /scratch/ref/gdlab/blast_db/nt_2022_03_08/nt \
    
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi