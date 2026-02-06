#!/bin/bash

#===============================================================================
# Name         : s06.1_bbmap.sh
# Description  : This script will run bbmap to estimate coverage for an assembly
# Usage        : sbatch s06.1_bbmap.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.2
# Created On   : 2020-04-17 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-09-13
#===============================================================================

#SBATCH --job-name=bbmap
#SBATCH --array=1-160%20
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm/bbmap/z_bbmap_%a.out
#SBATCH --error=slurm/bbmap/z_bbmap_%a.out

eval $( spack load --sh bbmap )

basedir="$PWD"
readsdir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools"
refdir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/MetaSpades"
outdir="${basedir}/d06_assembly_qc/bbmap"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

# Make output directory
mkdir -p ${outdir}/${sample}

set -x
time bbwrap.sh \
    in1=${readsdir}/${sample}_FW_clean.fastq.gz,${readsdir}/${sample}_singletons_CLEAN.fastq.gz \
    in2=${readsdir}/${sample}_RV_clean.fastq.gz,null \
    ref=${refdir}/${sample}/scaffolds.fasta \
    nodisk=t \
    covstats=${outdir}/${sample}/scaffold_covg_stats.txt \
    t=${SLURM_CPUS_PER_TASK} \
    append &> ${outdir}/${sample}/summ_covg_stats.txt
RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occurred in ${sample}!"
    exit $RC
fi