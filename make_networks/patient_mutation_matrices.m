%% Put patient mutation data in form for network analysis
%   Only contains proteins that are matric of interacting proteins
clear;clc;close all
load adjacency_interactions
load CosmicIncFus_samples_tcga

num_patients = size(samples_tcga,1);

patient_data = samples_tcga(:,1);

for i = 1:num_patients
    % Get indices of patient mutations, if that protein is in interacting matrix
    mut_protein_ids = [];
    mut_protein_names = [];
    mutations = samples_tcga{i,2};
    for j = 1:height(mutations)
        gene = mutations{j,1};
        pos = find(strcmp(gene,A_index));
        if ~isempty(pos)
            mut_protein_ids = [mut_protein_ids ;pos];
            mut_protein_names = [mut_protein_names; gene];
        end
    end
    patient_data{i,2} = mut_protein_ids;
    patient_data{i,3} = mut_protein_names;
end

save patient_mutation_data patient_data A_index