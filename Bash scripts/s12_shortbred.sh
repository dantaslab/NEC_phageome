#!/bin/bash
#===============================================================================
# File Name    : s13_shortbred.sh
# Description  : This script will profile short read sequencing data
# Usage        : sbatch s11_shortbred.sh
# Author       : Bejan Mahmud
# Version      : 1.1
# Created On   : 2020-05-31
# Last Modified: 2024-08-29
# Modified By  : Kailun Zhang
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=shortbred
#SBATCH --cpus-per-task=2
#SBATCH --array=74,383
#SBATCH --mem=1G
#SBATCH --output=slurm/shortbred/x_shortbred_%A_%a.out
#SBATCH --error=slurm/shortbred/y_shortbred_%A_%a.err

eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate /ref/gdlab/software/envs/shortbred

# Basedir
basedir="$PWD"
#indir="${basedir}/d03_cleanreads/final"
outdir="${basedir}/d13_shortbred1"

# Make output directory
mkdir -p ${outdir}

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Zot_DNA_sample5_list.txt`
indir=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240421_DNA_sample_path_list.txt`

# Start debug mode (will send commands to outfile)
set -x
#gunzip ${indir}/${sample}_.fastq.gz
#gunzip ${indir}/${sample}_R2.fastq.gz
time shortbred_quantify.py --markers /ref/gdlab/data/shortbred_markers/220529_shortbred_FMG-CARD-NCBIAMR_clustId-95pc.faa \
 --wgs ${indir}_R1_CLEAN.fastq.gz ${indir}_R2_CLEAN.fastq.gz \
 --results ${outdir}/${sample}_ShortBRED.txt \
 --tmp tmp_shortbred_reads/${sample}_quantify

rm -r tmp_shortbred_reads/${sample}_quantify

# Save error code for command
RC=$?

# End debug mode (will stop sending commands to outfile)
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occured for Sample ${sample}"
  exit $RC
fi