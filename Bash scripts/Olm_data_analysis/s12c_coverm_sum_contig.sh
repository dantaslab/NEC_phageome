#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm_sum.sh
# Description  : Summarize coverm data into a single file
# Usage        : sbatch s12c_coverm_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Modified On  : 2024-11-18
#===============================================================================

basedir="/scratch/gdlab/kailun/Phageome/NEC_Olm"
indir="${basedir}/d12_bowtie/CoverM_contig"
samplelist="${basedir}/d12_bowtie/sample_list.txt"
outfile="${indir}/coverm_sum_OlmSRR_95_contig.tsv"

# create output file
touch ${outfile}
#header="Bin\tGenome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Bases\tCovered Fraction"
header="Sample\tGenome\tMean\tTrimmed Mean\tCovered Bases\tVariance\tLength\tCount\tReads_per_base\tRPKM\tTPM"

echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        #assembly="${sample}"
        info=`grep -v "Contig" ${indir}/${sample}_coverm_output.tsv | awk -v name=${sample} '{ print name "\t" $0 }'`
        #echo -e "${assembly}\t${info}" >> ${outfile}
        echo -e "${info}" >> ${outfile}
done < ${samplelist}
