function [total_repression,feature_table] = predict(mRNA_ID,miRNA_seq,with_features)

%Calculating features and repression for all canonical sites for a
%specified miRNA in a specified mRNA

%mRNA_ID = The mRNA ENST ID [char]

%miRNA_seq = The full miRNA sequence [char]

%with_features - 0/1 index for including calculated features in the output or only the
%predicted repression [double]

addpath('Utils')
load('Data/gene_list.mat','gene_list')

gene_ind = find(strcmp(gene_list.ENST,mRNA_ID),1);

if isempty(gene_ind)
    error('Sorry, this ENST doesn''t exist in the database. predict_site may be used for custom genes')
end

load(sprintf('Data/genes_%d.mat',gene_list.ind(gene_ind)),'genes')

genes = genes(strcmp(genes.ENST,mRNA_ID),:);

UTR5 = genes.UTR5{1};
ORF = genes.ORF{1};
UTR3 = genes.UTR3{1};
RNA = [UTR5,ORF,UTR3];

s = genes.missing_len(1);
phastcons100{1} = [zeros(s,1); genes.phastcons100{1}];
phastcons20{1} = [zeros(s,1); genes.phastcons20{1}];
phylops100{1} = [zeros(s,1); genes.phylops100{1}];
phylops20{1} = [zeros(s,1); genes.phylops20{1}];
riboseq = genes.riboseq{1};

clear genes

ORF_start = length(UTR5)+1;
UTR3_start = length(UTR5)+length(ORF)+1;

site_table = find_potential_targets(RNA,miRNA_seq,ORF_start,UTR3_start,1);
site_table(site_table.RNA_start <= length(UTR5) + 2 * length(ORF) / 3,:) = [];

for i = 1:height(site_table)
    feature_row = predict_site(RNA,ORF_start,UTR3_start,miRNA_seq,site_table.RNA_start(i),...
        site_table.seed_type{i},phastcons20,phastcons100,phylops20,phylops100,riboseq,1);
    if i == 1
        feature_table = feature_row;
    else
        feature_table = [feature_table; feature_row];
    end
end

if ~with_features
    feature_table = feature_table(:,'repression');
end

feature_table = [site_table, feature_table];
total_repression = sum(feature_table.repression);