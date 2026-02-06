#!/bin/bash
#===============================================================================
#
# File Name    : s11b_busco_sum.sh
# Description  : Summarize busco data into a single file
# Usage        : sbatch s11b_busco_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 2.0
# Created On   : Nov 05 2021
# Modified On  : 2024-09-30
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=busco
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=slurm/x_busco_%A.out
#SBATCH --error=slurm/y_busco_%A.out

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm/BUSCO"
indir="${basedir}"
samplelist="${indir}/240928_Olm_srr_phage_pass3_list_ct3.txt"
outfile="${basedir}/240924_srr_busco_sum_ct3.tsv"

# create output file
touch ${outfile}
header="Assembly\t Complete (single)\t Complete (duplicated)\t Fragmented\t Missing"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        assembly="${sample}"
        single=`grep single-copy ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -f2`
        duplicate=`grep duplicated ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -f2`
        fregament=`grep Fragmented ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -f2`
        miss=`grep Missing ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -f2`
        echo -e "${assembly}\t${single}\t${duplicate}\t${fregament}\t${miss}" >> ${outfile}

done < ${samplelist}
