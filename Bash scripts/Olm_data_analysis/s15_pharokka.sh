#!/bin/bash
#===============================================================================
# File Name    : s08_pharokka.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s08_pharokka.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 2
# Modified     : 2022-07-20 
# Created      : 2024-06-27
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=pharokka
#SBATCH --array=1-5762%30
#SBATCH --mem=8G
#SBATCH --cpus-per-task=18
#SBATCH --output=slurm/x_pharokka_%a.out
#SBATCH --error=slurm/y_pharokka_%a.err

eval $( spack load --sh /buvzc6u )

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm/Pharokka/Pharokka_Olm_bac_phage_ct2_5762"
indir="${basedir}/single_contigs"
outdir="${basedir}"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240620_Olm_bac_phage_pass4_list_ct2.txt`

set -x

time 
python /ref/gdlab/software/envs/Mamba/Mambaforge-Linux-x86_64/envs/pharokka/bin/pharokka.py -i ${indir}/${sample}.fasta -o ${outdir}/${sample} -d /ref/gdlab/software/envs/Mamba/Mambaforge-Linux-x86_64/envs/pharokka -t ${SLURM_CPUS_PER_TASK}
python /ref/gdlab/software/envs/Mamba/Mambaforge-Linux-x86_64/envs/pharokka/bin/pharokka_plotter.py -i ${outdir}/${sample}/phanotate.ffn -n "${sample}" -o ${outdir}/${sample} --interval 8000 --annotations 0.5 --plot_title '${sample}' --truncate 50
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi