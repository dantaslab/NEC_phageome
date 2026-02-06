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
#SBATCH --array=101-160%10
#SBATCH --cpus-per-task=12
#SBATCH --mem=240G
#SBATCH --output=slurm/deconseq/x_deconseq_%a.out
#SBATCH --error=slurm/deconseq/y_deconseq_%a.err

#load the module
#module purge
#module load deconseq/0.4.3-chr38

eval $( spack load --sh deconseq-standalone ) 

#store the base directory
basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm" #$PWD
#store the input path   
indir="${basedir}/d00_Coassemblies"
#store the output path
outdir="${basedir}/d02_deconseq"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

#deconseq command

set -x
cp /lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/Coassemblies/${sample}_1.fastq.gz ${indir}
cp /lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/Coassemblies/${sample}_2.fastq.gz ${indir}

gunzip ${indir}/${sample}_1.fastq.gz
time deconseq-multi.pl \
     -f ${indir}/${sample}_1.fastq \
     -out_dir ${outdir} \
     -id ${sample}_fwd \
     -dbs hs_ref,cow\
     -t ${SLURM_CPUS_PER_TASK}
rm ${indir}/${sample}_1.fastq
RC_FW=$?

gunzip ${indir}/${sample}_2.fastq.gz
time deconseq-multi.pl \
     -f ${indir}/${sample}_2.fastq \
     -out_dir ${outdir} \
     -id ${sample}_rev \
     -dbs hs_ref,cow \
     -t ${SLURM_CPUS_PER_TASK}
rm ${indir}/${sample}_2.fastq
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