#!/bin/bash          
#===============================================================================      
# Name         : s07a_CenoteTaker3_prune.sh
# Description  : Discover and annoate virus sequences/genomes    
# Usage        : sbatch s07a_CenoteTaker3_prune.sh
# Author       : Kailun Zhang, kailun@wustl.edu         
# Version      : 1.0 
# Created On   : 2024-03-07
# Modified On  : 2024-09-08
#===============================================================================      
#Submission script for HTCF
#SBATCH --job-name=CT3 
#SBATCH --array=145-160%7
#SBATCH --mem=24G      
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm/x_ct3_%A_%a.out  
#SBATCH --error=slurm/y_ct3_%A_%a.err

eval $( spack load --sh /buvzc6u ) # miniconda3@4.10.3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
source activate /ref/gdlab/software/envs/Mamba/Mambaforge-Linux-x86_64/envs/ct3_env

basedir="${PWD}"
indir="${basedir}/d05_metaspades_coassembly"
#outdir="${basedir}/d10a_megahit_ct3"
#mkdir -p ${outdir}

# Read in the slurm array task                               
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`
runtitle=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`

set -x

cenotetaker3 -c ${indir}/${sample}.fasta -r ${runtitle} -p T -t ${SLURM_CPUS_PER_TASK}
#rm -r ${basedir}/d07_cenotetaker3/${runtitle}/sequin_and_genome_maps
#rm -r ${basedir}/d07_cenotetaker3/${runtitle}/ct_processing

RC=$?

set +x

#output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi