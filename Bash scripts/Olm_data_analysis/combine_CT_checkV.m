clear all;
close all;

[num, txt, Data] = xlsread('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Olm/Results_SRRs/240913_CheckV_quality_summary.xlsx');
checkV_info      = Data([2:size(txt,1)],:);
Contig_source    = txt([2:size(txt,1)],1);

Contig_query     = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Olm/Results_SRRs/240913_Olm_SRR_CT3sum.xlsx', 'Range','A2:A22457');

CT_query = [];
for i = 1:size(Contig_query,1)
    i
    if find(string(Contig_source)==string(Contig_query(i)))
        CT_query = [CT_query; checkV_info(string(Contig_source)==string(Contig_query(i)),:)];
    else
        Contig_query(i)
    end
end

writecell(CT_query,'240915_combine_CT3_checkV.txt');