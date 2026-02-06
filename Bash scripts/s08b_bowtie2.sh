#!/bin/bash
#===============================================================================
#
# File Name    : s12b_bowtie2.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12b_bowtie2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: May 31 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtie2
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-2120
#SBATCH --output=slurm/bowtie/x_bowtie2_%A_%a.out
#SBATCH --error=slurm/bowtie/y_bowtie2_%A_%a.out

module load bowtie2/2.3.4.1
module load samtools/1.6

basedir="/scratch/gdlab/kailun/Phageome/NEC"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/sample_list.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d12_bowtie/ref_list.txt`
indir="${basedir}/NEC_cleanreads_time"
#index is where the output from bowtie2-build is stored
index="${basedir}/d12_bowtie/d12_bowtieBuild/${ref}/${ref}"
outdir="${basedir}/d12_bowtie/d12_bowtie2/${sample}"

#make output directory
mkdir -p ${outdir}

time bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${index} \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 ${indir}/${sample}_R1_CLEAN.fastq.gz \
    -2 ${indir}/${sample}_R2_CLEAN.fastq.gz

# convert sam to bam, then generate sorted bam file and index it.
samtools view -bS ${outdir}/${sample}_mappedreads.sam > ${outdir}/${sample}_raw.bam
samtools sort ${outdir}/${sample}_raw.bam > ${outdir}/${sample}_sorted.bam
samtools index ${outdir}/${sample}_sorted.bam
# generate stats
samtools depth -aa ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_depth.txt
samtools stats ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_stats.txt

rm ${outdir}/${sample}_mappedreads.sam
rm ${outdir}/${sample}_raw.bam

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
