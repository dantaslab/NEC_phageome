#!/bin/bash

#===============================================================================
# File Name    : s13_metaphlan.sh
# Description  : Profiles taxonomic composition of a read set using MetaPhlAn
# Usage        : sbatch s13_metaphlan.sh
# Author       : Kailun Zhang
# Version      : 1.2
# Created On   : 2021-01-06 by Luke Diorio-Toth
# Last Modified: 2024-03-13
#===============================================================================

#SBATCH --job-name=metaphlan
#SBATCH --array=1-1158%25
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=slurm/metaphlan/x_metaphlan_sub_%a_%A.out
#SBATCH --error=slurm/metaphlan/y_metaphlan_sub_%a_%A.out

# There are multiple versions of metaphlan installed, which use different dbs
# run "spack find py-metaphlan" to see all of them
eval $( spack load --sh miniconda3@4.10.3 )
eval $( spack load --sh py-metaphlan@4.0.6 )

basedir="$PWD"
indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools/d03_bbtools_single"
outdir="${basedir}/d13_metaphlan"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

mkdir -p ${outdir}

set -x
time metaphlan \
  ${indir}/${sample}_FW_clean.fastq.gz,${indir}/${sample}_RV_clean.fastq.gz \
  --input_type fastq \
  --bowtie2out ${outdir}/${sample}.bowtie2.bz2 \
  -o ${outdir}/${sample}_profile.txt \
  --nproc ${SLURM_CPUS_PER_TASK} --bowtie2db /ref/gdlab/data/metaphlan4_db_Oct2022
  
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi