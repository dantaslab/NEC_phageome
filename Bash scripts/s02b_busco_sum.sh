#!/bin/bash
#===============================================================================
#
# File Name    : s01-2_busco_sum.sh
# Description  : Summarize busco data into a single file
# Usage        : sbatch s01-2_busco_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 05 2021
# Modified On  : May 26 2022
#===============================================================================

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/BUSCO_concate"
samplelist="${basedir}/VAMB/Concate_bins/concate_list.txt"
outfile="${indir}/BUSCO_summaries/busco_sum.tsv"

# create output file
touch ${outfile}
header="Assembly\t Complete (single)\t Complete (duplicated)\t Fragmented\t Missing"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        assembly="${sample}"
        single=`grep single-copy ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -c 2`
        duplicate=`grep duplicated ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -c 2`
        fregament=`grep Fragmented ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -c 2`
        miss=`grep Missing ${indir}/${sample}/short_summary.specific.bacteria_odb10.${sample}.txt | cut -c 2-4`
        echo -e "${assembly}\t${single}\t${duplicate}\t${fregament}\t${miss}" >> ${outfile}

done < ${samplelist}
