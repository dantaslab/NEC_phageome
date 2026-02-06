#!/bin/bash
#===============================================================================
#
# File Name    : s01_Cenote-Taker2_sum.sh
# Description  : Summarize busco data into a single file
# Usage        : s01_Cenote-Taker2_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 16 2021
# Modified On  : Nov 16 2021
#===============================================================================

basedir="/scratch/gdlab/kailun/Phageome/NEC"
indir="${basedir}/CenoteTaker2_nt"
CT2dir="${indir}"
samplelist="${basedir}/scaffold_list_new.txt"
outfile="${indir}/CT2_sum.tsv"

# create output file
touch ${outfile}

header="ORIGINAL_NAME\tCENOTE_NAME\tORGANISM_NAME\tEND_FEATURE\tLENGTH\tORF_CALLER\tNUM_HALLMARKS\tHALLMARK_NAMES\tBLASTN_INFO"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        CT_Data=`grep -P "${sample}" ${CT2dir}/${sample}/${sample}_CONTIG_SUMMARY.tsv`
        echo -e "${CT_Data}" >> ${outfile}

done < ${samplelist}
