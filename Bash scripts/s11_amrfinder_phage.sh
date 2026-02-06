#!/bin/bash
#===============================================================================
# File Name    : s07b_amrfinder_phage.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s07b_amrfinder_phage
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : Jun 25 2019 by Luke Diorio-Toth
# Last Modified: May 20 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=amrfinder
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/amrfinder/x_amrfinder_phage_%A.out
#SBATCH --error=slurm/amrfinder/y_amrfinder_phage_%A.out

#module load perl
#module load amrfinder

#eval $( spack load --sh amrfinder )
eval $( spack load --sh /4docezg )
eval $( spack load --sh perl )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}"
outdir="${basedir}/d07b_Phage_amrfinder"

#make output directory and read in the slurm array task
mkdir -p ${outdir}
#sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${indir}/230427_viral_unbin_list.txt`

set -x
time amrfinder --plus \
	-n ${indir}/Pass4_viral_concate_contig2.fna \
	-o ${outdir}/230520_phage_AMRFinder_out.tsv \
	--threads ${SLURM_CPUS_PER_TASK}
RC=$?
set +x


if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi