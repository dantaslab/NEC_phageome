#!/bin/bash
#===============================================================================
#
# File Name    : s13a_bowtieBuild_RNA_allphage.sh
# Description  : This script will index a ref for Burrows-Wheeler alignment
# Usage        : sbatch s13a_bowtieBuild_RNA_allphage.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2024-10-31
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtieBuild
#SBATCH --array=1-144
#SBATCH --mem=4G
#SBATCH --output=slurm/bowtie_allphage/B_bowtieBuild_%A_%a.out
#SBATCH --error=slurm/bowtie_allphage/B_bowtieBuild_%A_%a.out

eval $( spack load --sh bowtie2 )

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/Bowtie_All_phage_RNA/separate_contigs"
outdir="${basedir}/Bowtie_All_phage_RNA/BowtieBuild"
mkdir -p ${outdir}

ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_144_new.txt`
#sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_new.txt`

set -x

#time grep -A 1 ">${ref}\|C${ref}" ${indir}/Pass4_viral_concate_contig2.fna > ${basedir}/d12_bowtie/${ref}_contigs.fasta
#sed -i '/--/d' ${ref}_contigs.fasta
#make output directory
mkdir -p ${outdir}/${ref}
bowtie2-build ${indir}/${ref}_contigs.fasta ${outdir}/${ref}/${ref}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi