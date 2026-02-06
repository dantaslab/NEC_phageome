clear all;
close all;
       
[num, txt, hmmData] = xlsread('hmmsearch_concate_contigs.xlsx');
hmmsize = size(hmmData,1);
hmm = hmmData([1:hmmsize],:);

[num, txt, viral_contig] = xlsread('concate_list.xlsx');
vpf_hit= [];
for i=1:size(viral_contig,1)
    vpf_hit =[vpf_hit size(find(hmm==string(viral_contig{i})),1)];
end

writematrix(vpf_hit','vpf_hits.txt');