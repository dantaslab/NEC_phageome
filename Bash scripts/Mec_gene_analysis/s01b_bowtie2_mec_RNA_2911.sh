#!/bin/bash
#===============================================================================
#
# File Name    : s01b_bowtie2_mec_RNA.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s01b_bowtie2_mec_RNA.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2024-11-02
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtie2
#SBATCH --mem=8G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-2911%7
#SBATCH --output=slurm/bowtie/R_bowtie2_%A_%a.out
#SBATCH --error=slurm/bowtie/R_bowtie2_%A_%a.out

eval $( spack load --sh bowtie2 )
#load samtools/1.14
eval $( spack load --sh /6p5wlkk )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/250618_mec_RNA_mapping_sample_list_2911.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/250618_mec_RNA_mapping_ref_list_2911.txt`
indir=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/250618_mec_RNA_mapping_samplepath_list_2911.txt`
#index is where the output from bowtie2-build is stored
index="${basedir}/mec/BowtieBuild_single/${ref}/${ref}"
outdir="${basedir}/mec/Bowtie_mec_RNA_2911/BowtieMapping_single/${ref}_${sample}"

#make output directory
mkdir -p ${outdir}

time bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${index} \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 ${indir}_R1_CLEAN.fastq.gz \
    -2 ${indir}_R2_CLEAN.fastq.gz

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
