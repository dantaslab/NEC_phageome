#!/bin/bash
#===============================================================================
#
# File Name    : s01_cluster.sh
# Description  : Cluster base on ANI
# Usage        : sbatch s01_cluster.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 18 2021
# Modified On  : 2024-04-02
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=8
#SBATCH --job-name=blast
#SBATCH --mem=32G
#SBATCH --output=slurm/cluster/x_blast_%A.out
#SBATCH --error=slurm/cluster/y_blast_%A.out


#module load ncbi-blast/2.6.0+
#eval $( spack load --sh blast-plus )
eval $( spack load --sh /u7ssbm4 )

basedir="/scratch/gdlab/kailun/Phageome/NEC"

indir="${basedir}"
outdir="${basedir}/d01_cluster_phage_3038"

# create a temp directory with outdir to deposit all the scaffold files
mkdir -p ${outdir}

makeblastdb -in ${indir}/240402_NEC_viral_contig_after_drop_3038.fna -dbtype nucl -out ${outdir}/Phage_contig_db/Virus_nt_db
blastn -query ${indir}/240402_NEC_viral_contig_after_drop_3038.fna -db ${outdir}/Phage_contig_db/Virus_nt_db -outfmt '6 std qlen slen' -max_target_seqs 100000 -out ${outdir}/240402_blast_out.tsv -num_threads 8
python anicalc.py -i ${outdir}/240402_blast_out.tsv -o ${outdir}/240402_ani.tsv
#awk '$NF >= 85' ${outdir}/230501_binned_ani.tsv > ${outdir}/230501_binned_ani95_tcov85.tsv
python aniclust.py --fna ${indir}/240402_NEC_viral_contig_after_drop_3038.fna --ani ${outdir}/240402_ani.tsv --out ${outdir}/240402_cluster_phage_3038.tsv --min_ani 95 --min_tcov 85 --min_qcov 0