#!/bin/bash

#===============================================================================
# File Name    : s06.2_quast.sh
# Description  : This script will use quast to calculate qc stats on assemblies
# Usage        : sbatch s06.2_quast.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : 2020-02-04 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2024-09-13
#===============================================================================

#SBATCH --job-name=quast
#SBATCH --array=1-160%20
#SBATCH --mem=32G      
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm/quast/x_quast_%a.out
#SBATCH --error=slurm/quast/y_quast_%a.out

eval $( spack load --sh /xy2yz7j ) # py-quast@5.2.0

#store the base directory
basedir="$PWD"
indir="${basedir}/d05_metaspades_coassembly"
outdir="${basedir}/d06_assembly_qc/quast"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

mkdir -p ${outdir}/${sample}

set -x
time quast.py ${indir}/${sample}.fasta -l ${sample} -o ${outdir}/${sample}
RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occurred in ${sample}!"
    exit $RC
fi