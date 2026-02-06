#!/bin/bash
#===============================================================================
# Name         : s04.1_fastqc.sh
# Description  : Runs basic QC on raw and trimmed illumina reads
# Usage        : sbatch s04.1_fastQC.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.4
# Created On   : 2019_01_08 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-08-21
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=fastqc
#SBATCH --array=1-1158%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/fastqc/x_fastqc_%a_%A.out
#SBATCH --error=slurm/fastqc/y_fastqc_%a_%A.out

# load latest version of fastqc on spack (fastqc@0.11.9)
eval $( spack load --sh /ldosddd )

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
cleanin="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools/d03_bbtools_single"
cleanout="${basedir}/d04_read_qc_single/fastqc"
mkdir -p ${cleanout}
export JAVA_ARGS="-Xmx8000M"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

set -x

time fastqc ${cleanin}/${sample}_FW_clean.fastq.gz \
            ${cleanin}/${sample}_RV_clean.fastq.gz \
            -o ${cleanout} \
            -t ${SLURM_CPUS_PER_TASK}

RC=$?
set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
