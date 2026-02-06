#!/bin/bash
#===============================================================================
# File Name    : s06.4_assembly_qc.sh
# Description  : Consolidates output from s06.1 to s06.3
# Usage        : sbatch s06.4_assembly_qc.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2020-07-02 byLuke Diorio-Toth, ldiorio-toth@wustl.edu
# Modified On  : 2024-03-11
#===============================================================================

basedir="$PWD"
indir="${basedir}/d06_assembly_qc"
quastdir="${indir}/quast"
bbmapdir="${indir}/bbmap"
checkmdir="${indir}/checkm"
samplelist="${basedir}/240821_coassembly_list_160.txt"
outfile="${indir}/240915_assembly_qc.tsv"

# create output file
touch ${outfile}
header="Assembly\t# contigs (>= 0 bp)\t# contigs (>= 1000 bp)\t# contigs (>= 5000 bp)\t# contigs (>= 10000 bp)\t# contigs (>= 25000 bp)\t# contigs (>= 50000 bp)\tTotal length (>= 0 bp)\tTotal length (>= 1000 bp)\tTotal length (>= 5000 bp)\tTotal length (>= 10000 bp)\tTotal length (>= 25000 bp)\tTotal length (>= 50000 bp)\t# contigs\tLargest contig\tTotal length\tGC (%)\tN50\tN75\tL50\tL75\t# N's per 100 kbp\tAverage covg\tCovg SD\t% reads mapped\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
    # find assembly stats
    quast_stats=`find ${quastdir}/${sample} -name 'transposed_report.tsv' -exec sed '2!d' {} \;`
    covg=`grep "Average coverage" ${bbmapdir}/${sample}/summ_covg_stats.txt | cut -f 2`
    covg_sd=`grep "Standard deviation" ${bbmapdir}/${sample}/summ_covg_stats.txt | cut -f 2`
    pct_mapped=`grep "Percent mapped" ${bbmapdir}/${sample}/summ_covg_stats.txt | cut -f 2`
    checkm_stats=`grep -P "${sample}\t" ${checkmdir}/output.txt | cut -f 1 --complement`
    echo -e "${quast_stats}\t${covg}\t${covg_sd}\t${pct_mapped}\t${checkm_stats}" >> ${outfile}
done < ${samplelist}