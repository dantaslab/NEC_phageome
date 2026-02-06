#!/bin/bash
#===============================================================================
#
# File Name    : s02a_separate_contig.sh
# Description  : This script will separate contigs based on headers
# Usage        : sbatch s02a_separate_contig.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-05-21
# Last Modified: 2022-05-21
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=separate
#SBATCH --array=1-144
#SBATCH --mem=4G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/x_cp_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/y_cp_%A_%a.out

# awk 'BEGIN{while((getline<"None_end_contig_list.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f' all_seq.fna > None_end_contig.fna
# awk 'BEGIN{while((getline<"DTR_ITR_contig_list.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f' all_seq.fna > DTR_ITR_contig.fna

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/VAMB/all_ct2_contigs/None_end_contigs"
outdir="${indir}/separate"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_144_new.txt`

set -x

time grep -A 1 ">${sample}" ${indir}/None_end_contig.fna > ${outdir}/${sample}_none_end_contigs.fna

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi