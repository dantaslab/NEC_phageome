#!/bin/bash

#===============================================================================
# Name         : s01_trimmomatic.sh
# Description  : Trims adapter sequences and low-quality bases from illumina
#              : reads, producing a "cleaned" set of paried reads
# Usage        : sbatch s01_trimmomatic.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.3
# Created On   : 2019_01_08 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-03-11
#===============================================================================

#SBATCH --job-name=trimmomatic
#SBATCH --array=1-12
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --output=slurm/trimmomatic/x_trim_%a_%A.out
#SBATCH --error=slurm/trimmomatic/y_trim_%a_%A.out

eval $( spack load --sh trimmomatic )

basedir="/scratch/gdlab/kailun/Phageome/Pak_bulk"
indir="/lts/gdlab/seq5/240301_Pak_Bulk_pilot"
outdir="${basedir}/d01_clean_reads"

mkdir -p ${outdir}

# need to declare memory explicitely
export JAVA_ARGS="-Xmx1000M"

# choose which adapters to use:
# NexteraPE-PE.fa,
# TruSeq2-PE.fa, TruSeq2-SE.fa,
# TruSeq3-PE-2.fa, TruSeq3-PE.fa, TruSeq3-SE.fa
adapt="/ref/gdlab/data/trimmomatic_adapters/0.39/NexteraPE-PE.fa"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240311_sample_list.txt`

# R1 = fwd and R2 = rev
# P = paired, UP = unpaired
set -x
time trimmomatic \
    PE \
    -phred33 \
    -trimlog \
    ${outdir}/Paired_${sample}_trimlog.txt \
    ${indir}/${sample}_R1.fastq.gz \
    ${indir}/${sample}_R2.fastq.gz \
    ${outdir}/${sample}_FW_clean.fastq.gz \
    ${outdir}/${sample}_FW_clean_UP.fastq.gz \
    ${outdir}/${sample}_RV_clean.fastq.gz \
    ${outdir}/${sample}_RV_clean_UP.fastq.gz \
    ILLUMINACLIP:${adapt}:2:30:10:1:true \
    SLIDINGWINDOW:4:20 \
    LEADING:10 \
    TRAILING:10 \
    MINLEN:60
RC=$?
set +x

# combine unpaired reads into a single fastq and delete redundant files
cat ${outdir}/${sample}_FW_clean_UP.fastq.gz \
    ${outdir}/${sample}_RV_clean_UP.fastq.gz \
    > ${outdir}/${sample}_UP_clean.fastq.gz
rm ${outdir}/${sample}_FW_clean_UP.fastq.gz
rm ${outdir}/${sample}_RV_clean_UP.fastq.gz 

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
