#!/bin/bash
#===============================================================================
#
# File Name    : s11e_quast_sum.sh
# Description  : Consolidates output from s06 assembly qc tools
# Usage        : sbatch 11e_quast_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2024-03-08
# Modified On  : 2024-03-10
#===============================================================================

#SBATCH --job-name=quast
#SBATCH --mem=1G
#SBATCH --output=slurm/quast/z_quast_sum_%A.out
#SBATCH --error=slurm/quast/z_quast_sum_%A.out

basedir="$PWD"
indir="${basedir}/d11d_assembly_quast_metaspades_all"
quastdir="${indir}"
samplelist="${basedir}/240301_sample_list.txt"
outfile="${indir}/240310_assembly_quast_metaspades_all.tsv"

# create output file
touch ${outfile}
header="Assembly\t# contigs (>= 0 bp)\t# contigs (>= 1000 bp)\t# contigs (>= 5000 bp)\t# contigs (>= 10000 bp)\t# contigs (>= 25000 bp)\t# contigs (>= 50000 bp)\tTotal length (>= 0 bp)\tTotal length (>= 1000 bp)\tTotal length (>= 5000 bp)\tTotal length (>= 10000 bp)\tTotal length (>= 25000 bp)\tTotal length (>= 50000 bp)\t# contigs\tLargest contig\tTotal length\tGC (%)\tN50\tN90\tauN\tL50\tL90\t# N's per 100 kbp"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # find assembly stats
        quast_stats=`find ${quastdir}/${sample} -name 'transposed_report.tsv' -exec sed '2!d' {} \;`
        echo -e "${quast_stats}" >> ${outfile}
done < ${samplelist}