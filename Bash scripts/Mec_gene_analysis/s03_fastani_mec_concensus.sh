#!/bin/bash
#===============================================================================
# File Name    : s03_fastani_mec_concensus.sh
# Description  : Fast all vs. all ANI calculation of many genomes
# Usage        : sbatch s03_fastani_mec_concensus.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Modified     : 2024-10-24
# Created      : 2022-10-12 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=all-fastani
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm/fastani/z_fastani_%A.out
#SBATCH --error=slurm/fastani/z_fastani_%A.out

eval $( spack load --sh fastani )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/mec/mec_genes"
outdir="${basedir}/mec/fastANI"

# Make output directory
mkdir -p ${outdir}

# Find assemblies and write full path to a file
# for an all v all analysis, that file will be ref and query
find ${indir} -type f -name "*.fasta" >> ${outdir}/genome_paths.txt


set -x
time fastANI --ql ${outdir}/genome_paths.txt \
            --rl ${outdir}/genome_paths.txt \
            -o ${outdir}/fastani_out.txt \
            -t ${SLURM_CPUS_PER_TASK}  --fragLen 0 --minFraction 0
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occurred!"
  exit $RC
fi