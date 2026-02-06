clear all;
close all;
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);

[num, txt, Data] = xlsread('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx','ref');
Kingdom          = txt([2:size(txt,1)],2);
Phylum           = txt([2:size(txt,1)],3);
Class            = txt([2:size(txt,1)],4);
Order            = txt([2:size(txt,1)],5);
Family           = txt([2:size(txt,1)],6);
%Genus            = txt([2:size(txt,1)],7);
VC_source        = txt([2:size(txt,1)],11);

VC_status = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx', 'Sheet', 'sample', 'Range','I2:I4858');
VC_query  = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx', 'Sheet', 'sample', 'Range','J2:J4858');

Kingdom_query = [];
Phylum_query  = [];
Class_query   = [];
Order_query   = [];
Family_query  = [];
Genus_query   = [];
a             = [];
b             = [];
x = 0;
y = 0;
z = 0;
for i = 1:size(VC_query,1)
    if VC_status{i}=="Overlap"
        x = x+1;
        VC = regexp(VC_query{i}, "/", "split");
        if numel(unique(Family(ismember(VC_source, VC))))==1
            Kingdom_query = [Kingdom_query; unique(Kingdom(ismember(VC_source, VC)))];
            Phylum_query  = [Phylum_query; unique(Phylum(ismember(VC_source, VC)))];
            Class_query   = [Class_query; unique(Class(ismember(VC_source, VC)))];
            Order_query   = [Order_query; unique(Order(ismember(VC_source, VC)))];
            Family_query  = [Family_query; unique(Family(ismember(VC_source, VC)))];
%             Genus_query   = [Genus_query; Genus(find(string(VC_source)==VC{1}))];
        elseif numel(unique(Family(ismember(VC_source, VC))))==0
            Kingdom_query = [Kingdom_query;"Unassigned"];
            Phylum_query  = [Phylum_query; "Unassigned"];
            Class_query   = [Class_query; "Unassigned"];
            Order_query   = [Order_query; "Unassigned"];
            Family_query  = [Family_query; "Unassigned"];
        elseif numel(unique(Family(ismember(VC_source, VC))))>1
            Kingdom_query = [Kingdom_query;"Multi"];
            Phylum_query  = [Phylum_query; "Multi"];
            Class_query   = [Class_query; "Multi"];
            Order_query   = [Order_query; "Multi"];
            Family_query  = [Family_query; "Multi"];
        else
            a = [a VC_query(i)];
        end
    elseif VC_status{i}=="Clustered"
        y = y+1;
        if find(string(VC_source)==VC_query{i})
            Kingdom_query = [Kingdom_query; Kingdom(string(VC_source)==VC_query{i})];
            Phylum_query  = [Phylum_query; Phylum(string(VC_source)==VC_query{i})];
            Class_query   = [Class_query; Class(string(VC_source)==VC_query{i})];
            Order_query   = [Order_query; Order(string(VC_source)==VC_query{i})];
            Family_query  = [Family_query; Family(string(VC_source)==VC_query{i})];
%             Genus_query   = [Genus_query; Genus(find(string(VC_source)==string(VC_query(i))))];
        else
            Kingdom_query = [Kingdom_query;"Unassigned"];
            Phylum_query  = [Phylum_query; "Unassigned"];
            Class_query   = [Class_query; "Unassigned"];
            Order_query   = [Order_query; "Unassigned"];
            Family_query  = [Family_query; "Unassigned"];
        end
    elseif VC_status{i}=="Outlier" || VC_status{i}=="Singleton"
        z = z+1;
            Kingdom_query = [Kingdom_query;"Unassigned"];
            Phylum_query  = [Phylum_query; "Unassigned"];
            Class_query   = [Class_query; "Unassigned"];
            Order_query   = [Order_query; "Unassigned"];
            Family_query  = [Family_query; "Unassigned"];
%             Genus_query   = [Genus_query; string("Unassigned")];
    end
end

All_taxa = [Kingdom_query Phylum_query Class_query Order_query Family_query Genus_query];
% figure;
% pie(categorical(Family_query));
%writematrix(All_taxa,'Taxa_pass4_contig.txt');

%% Extract family previous analysis

% [num, txt, Data] = xlsread('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/busco_sum_0.05_contig.xlsx');
% name             = txt([2:size(txt,1)],23);
% phage_family     = txt([2:size(txt,1)],31);
% 
% sample_query = readcell('/Users/Kellen_1/Google
% Drive/WashU/Phageome/NEC/Data
% analysis/Concate_output/busco_sum_try1.xlsx', 'Sheet', 'Sheet1',
% 'Range','A1:A304'); % I've deleted this sheet; Basically paste the
% orgnism_name that is needed to be found
% 
% phage_family_query = [];
% for i = 1:size(sample_query,1)
%     if find(string(name)==string(sample_query(i)))
%         phage_family_query = [phage_family_query; phage_family(string(name)==string(sample_query(i)))];
%     else 
%         sample_query(i)
%     end
% end
% 
% % writecell(phage_family_query,'ct2_phage_family2.txt');

%% From Cenote-Taker2 results:
Contig_CT2    = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/busco_sum_try1.xlsx', 'Range','A2:A11095');
Family_CT2    = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/busco_sum_try1.xlsx', 'Range','Z2:Z11095');
Contig_target = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx','Sheet','sample' ,'Range','A2:A4858');
Family_target = [];
for j = 1:size(Contig_target,1)
    if find(string(Contig_CT2)==string(Contig_target(j)))
        Family_target = [Family_target; Family_CT2(string(Contig_CT2)==string(Contig_target(j)))];
    else
        Contig_target(j)
    end
end

% writecell(Family_target,'pass4_ct2_phagefamily.txt');

%% compile two method
family_vcontact2 = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx','Sheet','sample' ,'Range','O2:O4858');
family_ct2       = Family_target;

family_refined = [];
w = [];
v = [];
for x = 1:size(family_ct2,1)
    if family_vcontact2{x}=="Unassigned" || family_vcontact2{x}=="Multi"
        family_refined = [family_refined;string(family_ct2{x})];
    elseif string(family_vcontact2{x})==string(family_ct2{x})
        family_refined = [family_refined;string(family_ct2{x})];
    else
%         x 
%         Contig_target(x)
        family_refined = [family_refined;string(family_vcontact2{x})];
        w = [w string(family_vcontact2{x})];
        v = [v string(family_ct2{x})];
    end
end
% writematrix(family_refined,'pass4_refined_phagefamily.txt');

family_vcontact2 = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx','Sheet','sample' ,'Range','O2:O4858');
family_ct2       = Family_target;

family_refined = [];
w = [];
v = [];
for x = 1:size(family_ct2,1)
    if family_vcontact2{x}=="Unassigned" || family_vcontact2{x}=="Multi"
        family_refined = [family_refined;string(family_ct2{x})];
    elseif string(family_vcontact2{x})==string(family_ct2{x})
        family_refined = [family_refined;string(family_ct2{x})];
    else
%         x 
%         Contig_target(x)
        family_refined = [family_refined;string(family_vcontact2{x})];
        w = [w string(family_vcontact2{x})];
        v = [v string(family_ct2{x})];
    end
end
% writematrix(family_refined,'pass4_refined_phagefamily.txt');

family_refined(family_refined=="Unknown" | family_refined=="metagenomic plasmid" | family_refined=="Conjugative Transposon" | family_refined=="circular genetic element") = {'Others'};
% Family_combine = readcell('/Users/Kellen_1/Google Drive/WashU/Phageome/NEC/Data analysis/Concate_output/vContact2_genome_by_genome_refined.xlsx', 'Range','S2:S4858');
t = tiledlayout(1,3,'TileSpacing','compact');
ax1 = nexttile;
[n,y] = hist(categorical(family_vcontact2));
labels = arrayfun( @(ii) sprintf('%s (%.3f)', y{ii}, n(ii)./size(family_vcontact2,1)*100), 1:numel(n), 'uni', 0 );
pie(n, labels)
%pie(categorical(Family_query))
title('vConTACT','fontsize',20)

ax2 = nexttile;
[n,y] = hist(categorical(family_ct2));
labels = arrayfun( @(ii) sprintf('%s (%.3f)', y{ii}, n(ii)./size(family_ct2,1)*100), 1:numel(n), 'uni', 0 );
pie(n, labels)
%pie(categorical(Family_target))
title('Cenote-Taker2','fontsize',20)

ax3 = nexttile;
[n,y] = hist(categorical(family_refined));
labels = arrayfun( @(ii) sprintf('%s (%.2f)', y{ii}, n(ii)./size(family_refined,1)*100), 1:numel(n), 'uni', 0 );
pie(n, labels)
camroll(-90)
%pie(categorical(Family_combine))
title('Combine','fontsize',20)

