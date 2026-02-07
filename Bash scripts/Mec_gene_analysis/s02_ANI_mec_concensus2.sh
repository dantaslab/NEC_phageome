#!/bin/bash
#===============================================================================
#
# File Name    : s02_ANI_mec_concensus.sh
# Description  : Cluster base on ANI
# Usage        : sbatch s02_ANI_mec_concensus.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 18 2021
# Modified On  : 2025-04-21
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=2
#SBATCH --job-name=blast
#SBATCH --mem=2G
#SBATCH --output=slurm/cluster/x_blast_%A.out
#SBATCH --error=slurm/cluster/y_blast_%A.out


#module load ncbi-blast/2.6.0+
#eval $( spack load --sh blast-plus )
eval $( spack load --sh /u7ssbm4 )

basedir="/scratch/gdlab/kailun/Phageome/NEC"

indir="${basedir}/mec"
outdir="${basedir}/mec/ANI2"

# create a temp directory with outdir to deposit all the scaffold files
mkdir -p ${outdir}

makeblastdb -in ${indir}/250615_mec_consensus_refine.fna -dbtype nucl -out ${outdir}/mec_db/mec_nt_db
blastn -query ${indir}/250615_mec_consensus_refine.fna -db ${outdir}/mec_db/mec_nt_db -perc_identity 0 -evalue 1000 -outfmt 6 -out ${outdir}/250618_mec_blast_out.tsv -gapopen 5 -gapextend 2 -penalty -3 -reward 2
#python anicalc_ipd50.py -i ${outdir}/250618_mec_blast_out.tsv -o ${outdir}/250618_mec_blast_ani.tsv
#awk '$NF >= 85' ${outdir}/250421_ani.tsv > ${outdir}/250421_ani95_tcov85.tsv
#python aniclust.py --fna ${indir}/250421_NEC_9Imp_viral_contig_after_drop_early_407.fna --ani ${outdir}/250421_ani.tsv --out ${outdir}/250421_cluster_prune.tsv --min_ani 95 --min_tcov 85 --min_qcov 0