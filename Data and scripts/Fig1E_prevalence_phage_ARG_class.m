clear all;
close all;

[num_amr, txt_amr, Data_amr] = xlsread('Phage_pass4_AMRFinder_result.xlsx','AMR');
contig_AMR                   = txt_amr([2:size(txt_amr,1)],2);
AMR_class                    = txt_amr([2:size(txt_amr,1)],11);
AMR_class_list = readcell('Phage_pass4_AMRFinder_result.xlsx', 'Sheet', 'AMR_info', 'Range','A25:A35');

[num, txt, Data] = xlsread('NEC_phage_info.xlsx','Phage_all_samples');
%patient          = txt([2:size(txt,1)],1);
sample           = txt([2:size(txt,1)],3);
contig           = txt([2:size(txt,1)],4);

AMR_class_sample = [];
for j = 1:size(contig,1)
    if find(string(contig_AMR)==contig{j})
        idx_contig = find(string(contig_AMR)==contig{j});
        AMR_class_sample = [AMR_class_sample;strjoin(AMR_class(idx_contig),";")];
    else
        AMR_class_sample = [AMR_class_sample;"NA"];
    end
end

[num, txt, Data] = xlsread('NEC_metadata.xlsx');
sample_query     = txt([2:size(txt,1)],1);
date_after_birth = num(:,8);

% host_genus_all  = strjoin(host_genus,";");
% host_genus_list = regexp(host_genus_all, ";", "split");

info_all = [];
prevalence_all = [];
for i = 1:size(AMR_class_list,1)
    i
    taxa_idx  = contains(AMR_class_sample, AMR_class_list{i});
    sample_name = unique(sample(taxa_idx)); % unique sample names with this host genus
    PostAge = date_after_birth(contains(sample_query,sample_name));
    [n_pre date] = hist(PostAge,1:1:max(date_after_birth));
    [n_all date] = hist(date_after_birth,1:1:max(date_after_birth));
    prevalence = n_pre./n_all*100;
    info_all = [info_all; string(repelem(AMR_class_list(i),size(date,2))')];
    prevalence_all = [prevalence_all; [date' prevalence']];
    hold on;
end

writematrix(info_all,'Prevalence_phage_ARG_class.xlsx');
writematrix(prevalence_all,'Prevalence_phage_ARG_class.xlsx', 'Sheet',1,'Range','B1:C836')
