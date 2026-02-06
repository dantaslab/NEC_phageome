#!/bin/bash

#===============================================================================
# Name         : s04.2_multiqc.sh
# Description  : Consolidates output from fastqc.sh
# Usage        : sbatch s04.2_multiqc.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.4
# Created On   : 2019_01_08 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-08-21
#===============================================================================

#Submission script for HTCF
#SBATCH --job-name=multiqc
#SBATCH --mem=4G
#SBATCH --output=slurm/fastqc/z_multiqc_%A.out
#SBATCH --error=slurm/fastqc/z_multiqc_%A.out

eval $( spack load --sh py-multiqc )

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d04_read_qc_single"
cleanin="${indir}/fastqc"
cleanout="${indir}/multiqc"

set -x

time multiqc ${cleanin} -o ${cleanout} -n clean_multiqc

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi