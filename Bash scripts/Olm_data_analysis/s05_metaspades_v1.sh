#!/bin/bash
#===============================================================================
#
# File Name    : s05_metaspades.sh
# Description  : This script will assemble metagenomes with metaspades
# Usage        : sbatch s05_metaspades.sh
# Author       : Kailun Zhang
# Version      : 1.0
# Created On   : Wed Sep 11 15:58:08 CDT 2019 by Alaric D'Souza
# Last Modified: 2024-09-04
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=metaspades
#SBATCH --array=150
#SBATCH --cpus-per-task=20
#SBATCH --mem=640G
#SBATCH --output=slurm/metaspades/x_metaspades_%A_%a.out
#SBATCH --error=slurm/metaspades/y_metaspades_%A_%a.out

eval $( spack load --sh spades@3.15.3 )

#basedir if you want
basedir="$PWD"
indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools"
outdir="${basedir}/d05_metaspades_coassembly_killed2"

mkdir -p ${outdir}

ID=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

echo "Started: `date`"
echo "Host: `hostname`"

#start debug mode (will send commands to outfile)
set -x
#run and time command
time spades.py --memory ${SLURM_MEM_PER_NODE} --threads ${SLURM_CPUS_PER_TASK} --only-assembler --meta --tmp-dir /tmp/${ID} -1 ${indir}/${ID}_FW_clean.fastq.gz -2 ${indir}/${ID}_RV_clean.fastq.gz -o ${outdir}/${ID}
# --only-assembler 
#save error code for command
RC=$?
#exit debug mode
set +x

#output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi