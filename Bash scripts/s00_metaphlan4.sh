#!/bin/bash

#===============================================================================
# File Name    : s00_metaphlan.sh
# Description  : Profiles taxonomic composition of a read set using MetaPhlAn
# Usage        : sbatch s00_metaphlan.sh
# Author       : Kailun Zhang
# Version      : 1.2
# Created On   : 2021-01-06 by Luke Diorio-Toth
# Last Modified: 2024-04-23
#===============================================================================

#SBATCH --job-name=metaphlan
#SBATCH --array=1-2103%50
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=slurm/metaphlan/x_metaphlan_%A_%a.out
#SBATCH --error=slurm/metaphlan/y_metaphlan_%A_%a.out

# There are multiple versions of metaphlan installed, which use different dbs
# run "spack find py-metaphlan" to see all of them
eval $( spack load --sh py-metaphlan@4.0.6 )

basedir="$PWD"
indir=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240421_DNA_sample_path_list.txt`
outdir="${basedir}/Metaphlan4"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Zot_DNA_sample5_list.txt`

set -x
time metaphlan \
  ${indir}_R1_CLEAN.fastq.gz,${indir}_R2_CLEAN.fastq.gz \
  --input_type fastq \
  --bowtie2out ${outdir}/${sample}.bowtie2.bz2 \
  -o ${outdir}/${sample}_profile.txt \
  --nproc ${SLURM_CPUS_PER_TASK} \
  --bowtie2db /ref/gdlab/data/metaphlan4_db_Jun2023 #--unclassified_estimation
  
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
