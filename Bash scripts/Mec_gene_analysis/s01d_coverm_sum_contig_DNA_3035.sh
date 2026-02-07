#!/bin/bash
#===============================================================================
#
# File Name    : s01c_coverm_sum_contig_DNA.sh
# Description  : Summarize coverm data into a single file
# Usage        : sbatch s01c_coverm_sum_contig_DNA.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Modified On  : 2024-11-03
#===============================================================================

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/mec/Bowtie_mec_DNA_3035/CoverM_contig_single"
samplelist="${basedir}/250618_mec_DNA_mapping_ref_sample_list_3035.txt"
outfile="${indir}/250618_coverm_sum_mec_DNA_95_contig_3035.tsv"

# create output file
touch ${outfile}
#header="Bin\tGenome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Bases\tCovered Fraction"
header="Sample\tGenome\tMean\tTrimmed Mean\tCovered Bases\tVariance\tLength\tCount\tReads_per_base\tRPKM\tTPM"

echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        #assembly="${sample}"
        info=`grep -v "Contig" ${indir}/${sample}_coverm_DNA_output.tsv | awk -v name=${sample} '{ print name "\t" $0 }'`
        #echo -e "${assembly}\t${info}" >> ${outfile}
        echo -e "${info}" >> ${outfile}
done < ${samplelist}
