#!/bin/bash
#===============================================================================
# File Name    : s02b_vamb.sh
# Description  : Discover and annoate virus sequences/genomes
# Usage        : sbatch s02b_vamb.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.1
# Modified     : 2022-05-11
# Created      : 2022-05-13
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=VAMB_s
#SBATCH --cpus-per-task=8
#SBATCH --array=1-144
#SBATCH --mem=6G
#SBATCH --output=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/x_vamb_s_%A_%a.out
#SBATCH --error=/scratch/gdlab/kailun/Phageome/NEC/slurm/vamb/y_vamb_s_%A_%a.out

#module use /opt/htcf/modules
#module use /opt/htcf/modules-legacy
#module use /opt/apps/labs/gdlab/modules

module load miniconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Install vamb
# pip install https://github.com/RasmussenLab/vamb/archive/v3.0.3.zip

module load minimap2
module load samtools
#conda activate /opt/apps/labs/gdlab/software/samtools

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/NEC_cleanreads"
outdir="${basedir}/VAMB/bam_saparate4"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_144.txt`
sample_new=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_144_new.txt`

set -x
module load miniconda3
python concatenate.py ${basedir}/VAMB/catalogues/catalogue_${sample_new}.fna.gz ${basedir}/VAMB/all_ct2_contigs/final_combined_virus_sequences_${sample_new}.fna
minimap2 -d ${basedir}/VAMB/catalogues/catalogue_${sample_new}.mmi ${basedir}/VAMB/catalogues/catalogue_${sample_new}.fna.gz
minimap2 -t 8 -N 5 -ax sr ${basedir}/VAMB/catalogues/catalogue_${sample_new}.mmi ${indir}/${sample}_R1.fastq.gz ${indir}/${sample}_R2.fastq.gz | samtools view -F 3584 -b --threads 8 > ${outdir}/${sample_new}.bam

#module load miniconda3
python /home/kailun/.local/bin/vamb --outdir /scratch/gdlab/kailun/Phageome/NEC/VAMB/vamb_${sample_new} --fasta ${basedir}/VAMB/catalogues/catalogue_${sample_new}.fna.gz --bamfiles ${outdir}/${sample_new}.bam -t 2 -o C --minfasta 1


RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
