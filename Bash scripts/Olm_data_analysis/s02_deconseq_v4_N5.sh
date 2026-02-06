#!/bin/bash

#===============================================================================
# File Name    : s02_deconseq.sh
# Description  : This script will remove human and cow read contamination from fastq reads in parallel
# Usage        : sbatch s02_deconseq.sh
# Author       : Kailun Zhang
# Version      : 1.0
# Created On   : 2020-04-14 by Rhiannon Vargas
# Last Modified: 2024-08-21
#===============================================================================
# Submission script for HTCF

#SBATCH --job-name=deconseq
#SBATCH --array=601-869%50
#SBATCH --cpus-per-task=24
#SBATCH --mem=72G
#SBATCH --output=slurm/deconseq/x_deconseq_N5_%a.out
#SBATCH --error=slurm/deconseq/y_deconseq_N5_%a.err

#load the module
#module purge
#module load deconseq/0.4.3-chr38

eval $( spack load --sh deconseq-standalone ) 
eval $( spack load --sh pigz )

#store the base directory
basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm" #$PWD
#store the input path   
indir="${basedir}/d01_clean_reads_single"
#store the output path
outdir="${basedir}/d02_deconseq_N5"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

#deconseq command

set -x
cp /lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/Rename/PRJNA417343/${sample}_1.fastq.gz ${indir}
cp /lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/Rename/PRJNA417343/${sample}_2.fastq.gz ${indir}

gunzip ${indir}/${sample}_1.fastq.gz
time deconseq-multi.pl \
     -f ${indir}/${sample}_1.fastq \
     -out_dir ${outdir} \
     -id ${sample}_fwd \
     -dbs hs_ref,cow\
     -t ${SLURM_CPUS_PER_TASK}
rm ${indir}/${sample}_1.fastq
pigz -p24 ${outdir}/${sample}_fwd_clean.fq
pigz -p24 ${outdir}/${sample}_fwd_cont.fq
RC_FW=$?

gunzip ${indir}/${sample}_2.fastq.gz
time deconseq-multi.pl \
     -f ${indir}/${sample}_2.fastq \
     -out_dir ${outdir} \
     -id ${sample}_rev \
     -dbs hs_ref,cow \
     -t ${SLURM_CPUS_PER_TASK}
rm ${indir}/${sample}_2.fastq
pigz -p24 ${outdir}/${sample}_rev_clean.fq
pigz -p24 ${outdir}/${sample}_rev_cont.fq
RC_RV=$?

set +x

RC=$((RC_FW + RC_RV))

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured!"
    exit $RC
fi