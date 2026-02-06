#!/bin/bash
#===============================================================================
#
# File Name    : s07c_CenoteTaker3_sum.sh
# Description  : Summarize CT2 CONTIG_SUMMARY.tsv into a single file
# Usage        : sbatch s07c_CenoteTaker3_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 16 2021
# Modified On  : 2024-09-13
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=sum_ct3
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=slurm/z_ct3sum_%A.out
#SBATCH --error=slurm/z_ct3sum_%A.out

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d07_cenotetaker3"
CT3dir="${indir}"
samplelist="${basedir}/240821_coassembly_list_160.txt"
outfile="${indir}/240913_Olm_SRR_CT3sum.tsv"

# create output file
touch ${outfile}

header="contig\tinput_name\torganism\tvirus_seq_length\tend_feature\tgene_count\tvirion_hallmark_count\trep_hallmark_count\tRDRP_hallmark_count\tvirion_hallmark_genes\trep_hallmark_genes\tRDRP_hallmark_genes\ttaxonomy_hierarchy\tORF_caller"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        CT_Data=`grep -P "${sample}" ${CT3dir}/${sample}/${sample}_virus_summary.tsv`
        echo -e "${CT_Data}" >> ${outfile}

done < ${samplelist}
