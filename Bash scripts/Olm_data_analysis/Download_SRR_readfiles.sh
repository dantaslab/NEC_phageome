eval $( spcack load --sh sratoolkit )
vdb-config -i


eval $( spack load --sh entrezdirect )
eval $( spack load --sh py-parallel-fastq-dump )
parallel-fastq-dump --tmpdir /tmp -O /scratch/gdlab/kailun/Phageome/NEC_Olm/SRA_Reads --threads 4 --split-files --gzip --sra-id SRR6262353  SRR6262316


29  SRR6257420 Have problem
SRR5963258
SRR5963000
SRR5420279

vdb-dump SRR6262351 -R1 -C READ_TYPE
vdb-dump SRR5963000 -R1 -C READ_TYPE


awk 'BEGIN{while((getline<"240620_Olm_bac_phage_pass3_list_ct2.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f' 240620_Olm_bac_CT2.fasta > 240620_Olm_bac_phage_pass3_ct2.fasta
eval $( spack load --sh py-pyfasta )
pyfasta split --header "%(seqid)s.fasta" 240620_Olm_bac_phage_pass3_ct2.fna

eval $( spack load --sh entrezdirect )
esearch -db sra -query PRJNA417343 \
       | efetch --format runinfo \
       | cut -d ',' -f 1 \
       | grep SRR > PRJNA417343_SRR_list.txt


