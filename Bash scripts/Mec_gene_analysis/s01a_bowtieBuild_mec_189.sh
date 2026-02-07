#!/bin/bash
#===============================================================================
#
# File Name    : s01a_bowtieBuild_mec.sh
# Description  : This script will index a ref for Burrows-Wheeler alignment
# Usage        : sbatch s01a_bowtieBuild_mec.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2024-10-31
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtieBuild
#SBATCH --mem=4G
#SBATCH --array=1-189
#SBATCH --output=slurm/bowtie/B_bowtieBuild_%A_%a.out
#SBATCH --error=slurm/bowtie/B_bowtieBuild_%A_%a.out

eval $( spack load --sh bowtie2 )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/mec/mec_genes/single"
outdir="${basedir}/mec/BowtieBuild_single"

ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/250615_phage_mec_protein_189.txt`

set -x

#time grep -A 1 ">${ref}\|C${ref}" ${indir}/Pass4_viral_concate_contig2.fna > ${basedir}/d12_bowtie/${ref}_contigs.fasta
#sed -i '/--/d' ${ref}_contigs.fasta
#make output directory
mkdir -p ${outdir}/${ref}
bowtie2-build ${indir}/${ref}.fasta ${outdir}/${ref}/${ref}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi