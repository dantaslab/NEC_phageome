#!/bin/sh
#===============================================================================
#
# File Name    : s16b_humann3_renorm.sh
# Description  : Normalize the abundance output files of HUMAnN 
# Usage        : sbatch s16b_humann3_renorm.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Modified     : Kailun Zhang
# Version      : 1.1
# Created On   : Mon Jan  02 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=humann_relab
#SBATCH --mem=2G
#SBATCH --array=1-1158%500
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm/humann/r_humann3-norm_%a.out
#SBATCH --error=slurm/humann/r_humann3-norm_%a.out
# load python and activate humann venv
#. /ref/gdlab/software/envs/humann3/bin/activate
source activate /ref/gdlab/software/envs/humann3

basedir="$PWD"
indir="${basedir}/d16_humann3"
outdir="${indir}/relab"
#make output directories
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

set -x

# Run Humann3 humann_renorm_table utility script
humann_renorm_table --input ${indir}/${sample}_genefamilies.tsv \
--units relab \
--output ${outdir}/${sample}_genefamilies_relab.tsv 

humann_renorm_table --input ${indir}/${sample}_pathabundance.tsv \
--units relab \
--output ${outdir}/${sample}_pathabundance_relab.tsv 

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi

