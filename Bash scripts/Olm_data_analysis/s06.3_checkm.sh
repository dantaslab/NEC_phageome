#!/bin/bash

#===============================================================================
# Name         : s06.3_checkm.sh
# Description  : This script will run checkm on a unicycler output dir
# Usage        : sbatch s06.3_checkm.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.6
# Created On   : 2020-08-22 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: 2024-09-13
#===============================================================================

#SBATCH --job-name=checkm
#SBATCH --cpus-per-task=16
#SBATCH --mem=320G
#SBATCH --output=slurm/checkm/x_checkm_%A.out
#SBATCH --error=slurm/checkm/y_checkm_%A.out

eval $( spack load --sh py-checkm-genome )

basedir="$PWD"
indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/MetaSpades"
tmpdir="${basedir}/tmp"
outdir="${basedir}/d06_assembly_qc/checkm"

mkdir -p ${outdir}

# copy scaffold files from unicycler to a temporary directory
# (checkM requires all your files being in a single folder)
# IMPORTANT: the find command is meant to find the output dirs from unicycler.
mkdir -p ${tmpdir}
for dir in `find ${indir} -type d -maxdepth 1 -mindepth 1`; do
  cp $dir/scaffolds.fasta $tmpdir/$(basename $dir).fasta
done

set -x
time checkm lineage_wf -f ${outdir}/output.txt \
    -t ${SLURM_CPUS_PER_TASK} \
    -x fasta \
    --tab_table \
    ${tmpdir} ${outdir}
RC=$?
set +x

# remove temp files
rm -r ${tmpdir}

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi