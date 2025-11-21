clear all;
close all;

[num, txt, Data] = xlsread('NEC_phage_info.xlsx','Phage_all_samples');
patient          = txt([2:size(txt,1)],1);
sample           = txt([2:size(txt,1)],3);
contig           = txt([2:size(txt,1)],4);
metabolism       = txt([2:size(txt,1)],5);

[num, txt, Data] = xlsread('NEC_metadata.xlsx');
sample_query     = txt([2:size(txt,1)],1);
date_after_birth = num(:,8);

host_metabolism_all  = strjoin(metabolism,";");
host_metabolism_list = regexp(host_metabolism_all, ";", "split");
host_metabolism_name = unique(host_metabolism_list);
% [string(unique(host_metabolism_list))' hist(categorical(host_metabolism_list),unique(host_metabolism_list))']
host_metabolism_name = host_metabolism_name(host_metabolism_name~="NA");

info_all = [];
prevalence_all = [];
for i = 1:size(host_metabolism_name,2)
    i
    taxa_idx  = contains(metabolism, host_metabolism_name{i});
    sample_name = unique(sample(taxa_idx)); % unique sample names with this metabolism
    PostAge = date_after_birth(contains(sample_query,sample_name));
    [n_pre date] = hist(PostAge,1:1:max(date_after_birth));
    [n_all date] = hist(date_after_birth,1:1:max(date_after_birth));
    prevalence = n_pre./n_all*100;
    info_all = [info_all; string(repelem(host_metabolism_name(i),size(date,2))')];
    prevalence_all = [prevalence_all; [date' prevalence']];
    hold on;
end

writematrix(info_all,'Prevalence_metabolism.xlsx');
writematrix(prevalence_all,'Prevalence_metabolism.xlsx', 'Sheet',1,'Range','B1:C912')

