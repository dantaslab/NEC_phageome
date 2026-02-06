#!/bin/bash

#===============================================================================
# Name         : s03_bbtools.sh
# Description  : This script will run bbtools repair to repair deconseq reads
# Usage        : sbatch s03_bbtools.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2019_01_21 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-10-01
#===============================================================================

#Submission script for HTCF
#SBATCH --job-name=bbtools
#SBATCH --array=148-297
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --output=slurm/bbtools/x_bbtools_%a.out
#SBATCH --error=slurm/bbtools/y_bbtools_%a.err

eval $( spack load --sh bbmap@38.63 )
eval $( spack load --sh pigz )


basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d02_deconseq_N2"
outdir="${basedir}/d03_bbtools_single"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

set -x

repair.sh ow=t in=${indir}/${sample}_fwd_clean.fq.gz \
                in2=${indir}/${sample}_rev_clean.fq.gz \
                out=${outdir}/${sample}_FW_clean.fastq \
                out2=${outdir}/${sample}_RV_clean.fastq \
                outs=${outdir}/${sample}_singletons_CLEAN.fastq \
                repair=t tossbrokenreads=t 

#gzip -f ${outdir}/${sample}_FW_clean.fastq
#gzip -f ${outdir}/${sample}_RV_clean.fastq

pigz -p20 ${outdir}/${sample}_FW_clean.fastq
pigz -p20 ${outdir}/${sample}_RV_clean.fastq
pigz -p20 ${outdir}/${sample}_singletons_CLEAN.fastq 

RC=$?

set +x

#output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi