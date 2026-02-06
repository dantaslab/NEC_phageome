#!/bin/sh
#===============================================================================
#
# File Name    : s16_humann3.sh
# Description  : Assigns taxonomic function to filtered and trimmed reads 
# Usage        : sbatch s16_humann3.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Modified     : Kailun Zhang on 2-24-11-27
# Version      : 1.1
# Created On   : Mon Jan  02 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=humann3
#SBATCH --mem=64G
#SBATCH --array=1-1158%25
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/humann/x_humann3_%a.out
#SBATCH --error=slurm/humann/y_humann3_%a.out

#eval $( spack load --sh miniconda3@4.10.3 )
#eval $( spack load --sh py-metaphlan@4.0.6 )

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh 
source activate /ref/gdlab/software/envs/humann3
eval $( spack load --sh glpk )

basedir="$PWD"
indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools/d03_bbtools_single"
outdir="${basedir}/d16_humann3"

#make output directories
mkdir -p ${outdir}

cat_outdir="${outdir}/concat"
mkdir -p ${cat_outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/241001_sample_list_1158.txt`

set -x

# Concatenate paired reads 
cat ${indir}/${sample}_FW_clean.fastq.gz ${indir}/${sample}_RV_clean.fastq.gz > ${cat_outdir}/${sample}_concat.fastq.gz

# Run Humann3 on concatenated reads fastq files
#humann --remove-temp-output --input ${cat_outdir}/${sample}_concat.fastq.gz --output ${outdir} --threads 8 --input-format fastq.gz --nucleotide-database /ref/gdlab/data/humann3_db/chocophlan --protein-database /ref/gdlab/data/humann3_db/uniref --taxonomic-profile ${basedir}/d13_metaphlan/${sample}_profile.txt

time humann \
        --input ${cat_outdir}/${sample}_concat.fastq.gz \
        --input-format fastq.gz \
        --output ${outdir} \
        --output-basename ${sample} \
        --threads ${SLURM_CPUS_PER_TASK} \
        --nucleotide-database /ref/gdlab/data/humann3_db/chocophlan \
        --protein-database /ref/gdlab/data/humann3_db/uniref \
        --taxonomic-profile ${basedir}/d13_metaphlan/${sample}_profile.txt \
        --remove-temp-output

rm ${cat_outdir}/${sample}_concat.fastq.gz

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi

