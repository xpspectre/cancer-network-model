%% Import and get only breast cancer data
% data = CosmicMutantExportCensusv68210114(strcmp('breast',CosmicMutantExportCensusv68210114.Primarysite),:);
% data = CosmicMutantExportIncFusv68210114(strcmp('breast',CosmicMutantExportIncFusv68210114.Primarysite),:);
% 
% save CosmicIncFus_breast data

load CosmicIncFus_breast

%% List of mutated genes
genes = unique(data.Genename);

%% List of sample names
samples = unique(data.Samplename);

%% Make a table for each sample that contains all the mutations for that sample
for i = 1:length(samples)
    samples{i,2} = table;
end

for i = 1:height(data)
    entry = data(i,:);
    name = entry.Samplename;
    sample_idx = find(strcmp(name,samples(:,1)));
    samples{sample_idx,2} = [samples{sample_idx,2}; entry];
end
save CosmicIncFus_samples samples

%% Get samples that have 2 or more entries/mutations
samples_multi = cell(0,2);
for i = 1:length(samples)
    if height(samples{i,2}) > 1
        samples_multi = [samples_multi;samples(i,:)];
    end
end
save CosmicIncFus_samples_multi samples_multi

%% Get TCGA samples only - these should be better quality/complete mutation sets
samples_tcga = cell(0,2);
for i = 1:length(samples)
    name = samples{i,1};
    if length(name) > 3 && strcmp(name(1:4),'TCGA')
        samples_tcga = [samples_tcga;samples(i,:)];
    end
end
save CosmicIncFus_samples_tcga samples_tcga

%% Get statistics of TCGA breast cancer data
counts = zeros(size(samples_tcga,1),1);
for i = 1:length(samples_tcga)
    counts(i) = height(samples_tcga{i,2});
end
mean_counts = mean(counts);

%% Get list of all genes mutated in TCGA breast cancer data
mutated_genes = {};
mutated_gene_ids = [];
mutated_gene_id_pairs = {};
for i = 1:length(samples_tcga)
    entries = samples_tcga{i,2};
    for j = 1:height(entries)
        entry = entries(j,:);
        gene_name = entry{1,1}{1};
        gene_id = entry{1,2};
        if ~ismember(gene_name,mutated_genes)
            mutated_genes = [mutated_genes; gene_name];
        end
        if ~isnan(gene_id) && ~ismember(gene_id,mutated_gene_ids)
            mutated_gene_ids = [mutated_gene_ids; gene_id];
            mutated_gene_id_pairs = [mutated_gene_id_pairs; {gene_id, gene_name}];
        end
    end
end
mutated_genes = sort(mutated_genes); % alphabetic order
mutated_gene_ids = sort(mutated_gene_ids);
mutated_gene_id_pairs = sortrows(mutated_gene_id_pairs,1);

save gene_id_list mutated_gene_id_pairs

T = cell2table(mutated_gene_id_pairs,'VariableNames',{'HGNC_ID','Name'});
writetable(T,'gene_id_list.txt','Delimiter','\t');