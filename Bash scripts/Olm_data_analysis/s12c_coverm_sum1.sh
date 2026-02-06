#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm_sum.sh
# Description  : Summarize coverm data into a single file
# Usage        : sbatch s12c_coverm_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Modified On  : 2022-06-04
#===============================================================================

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/d12_bowtie_unbin/d12_coverm_2"
samplelist="${basedir}/d12_bowtie_unbin/sample_list.txt"
outfile="${indir}/coverm_sum.tsv"

# create output file
touch ${outfile}
header="Bin\tGenome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Bases\tCovered Fraction\tVariance\tLength\tCount\tReads_per_base"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        assembly="${sample}"
        info=`grep contigs ${indir}/${sample}_coverm_output.tsv`
        echo -e "${assembly}\t${info}" >> ${outfile}
done < ${samplelist}
