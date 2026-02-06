clear all;
close all;
       
[num, txt, hmmData] = xlsread('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Olm/Results_SRRs/240915_hmmsearch_SRR_ct3_contig_list.xlsx');
hmmsize = size(hmmData,1);
hmm = hmmData([1:hmmsize],:);

% [num, txt, viral_contig] = xlsread('contig_list.xlsx');
[num, txt, Data] = xlsread('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Olm/Results_SRRs/240913_Olm_SRR_CT3_all_info.xlsx');
viral_contig = Data([2:size(Data,1)],1);
vpf_hit= [];
for i=1:size(viral_contig,1)
    i
    vpf_hit =[vpf_hit size(find(hmm==string(viral_contig{i})),1)];
end

writematrix(vpf_hit','240915_vpf_hits_ct3.txt');