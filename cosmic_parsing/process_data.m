%% Import and get only breast cancer data
% data = CosmicMutantExportCensusv68210114(strcmp('breast',CosmicMutantExportCensusv68210114.Primarysite),:);
data = CosmicMutantExportIncFusv68210114(strcmp('breast',CosmicMutantExportIncFusv68210114.Primarysite),:);

save CosmicIncFus_breast data

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
