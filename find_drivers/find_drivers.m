%% Get networks from hotnet
clear;clc;close all

f = fopen('components_multi.txt','r');

networks = {};
tline = fgetl(f);
while ischar(tline)
    C = strsplit(tline,'\t');
    networks = [networks;{C}];
    tline = fgetl(f);
end
fclose(f);

% Make a list of all genes in all networks
all_mutations = {};
for i = 1:size(networks,1)
    net_genes = networks{i};
    for j = 1:size(net_genes,2)
        all_mutations = [all_mutations; net_genes{j}];
    end
end

save all_mutations all_mutations

f = fopen('all_mutations.txt','w');
for i = 1:size(all_mutations,1)
    fprintf(f,'%s\n',all_mutations{i,1});
end
fclose(f);

%% Get mutations for patients
load patient_data

num_subnets = size(networks,1);
num_patients = size(patient_data,1);

%% Get matrices that hold which patient has which mutations for each subnet
net_matrices = cell(num_subnets,1); % holds gene/patient mutation data for each subnetwork
for i = 1:num_subnets
    network = networks{i,1};
    num_sub_genes = size(network,1);
    
    % Preallocate matrix of genes x patients with element = 1 if patient k
    % has mutation in gene j
    mutations = zeros(num_sub_genes,num_patients);
    
    for j = 1:num_patients
        p_genes = patient_data{j,2};
        for k = 1:length(p_genes)
            p_gene = p_genes{k};
            pos = find(strcmp(p_gene,network));
            if ~isempty(pos)
                mutations(pos,j) = 1;
            end
        end
    end
    net_matrices{i} = mutations;
end

%% Calculate mutual exclusivity
%     # x1 and x2 are vectors of 1 and 0, 1 for
%     # mutated and 0 for not mutated.
pairwise_mex = [];
for i = 1:num_subnets
    net_mat = net_matrices{i};
    num_sub_genes = size(net_mat,1);
    gene_counts = sum(net_mat,2); % Nx1 vector of number of times each gene in subnet is found in patients
    gene_probs = gene_counts/num_patients; % Nx1 vector of probabilities that each gene is found in patients
    
    % Get pairwise mutual exclusivities
    pair_positions = combnk(1:1:num_sub_genes,2);
    for j = 1:size(pair_positions,1)
        p1 = pair_positions(j,1);
        p2 = pair_positions(j,2);
        
        % Calculate mutual exclusivity
        joint_count = sum(net_mat(p1,:).*net_mat(p2,:)); % number of times mutation appears in both patients at once
        joint_prob = joint_count/num_patients;
        chi = joint_prob/(gene_probs(p1)*gene_probs(p2));
        
        % Calculate p-val
        pval = hygecdf(joint_count,num_patients,gene_counts(p1),gene_counts(p2));
        
        pairwise_mex = [pairwise_mex; i, p1, p2, chi, pval];
    end
end

%% Process results
% Sort by pval
pairwise_mex = sortrows(pairwise_mex,5);

labels = {'network','gene1','gene2','chi','pval'};

save pairwise_mex pairwise_mex labels networks

f = fopen('pairwise_mex.txt','w');
fprintf(f,'network\tgene1\tgene2\tchi\tpval\n');
for i = 1:size(pairwise_mex,1)
    % Get gene names
    network_genes = networks{pairwise_mex(i,1)};
    gene1 = network_genes{pairwise_mex(i,2)};
    gene2 = network_genes{pairwise_mex(i,3)};
    
    % Print to file
    fprintf(f,'%d\t%s\t%s\t%f\t%f\n',pairwise_mex(i,1),gene1,gene2,...
                                     pairwise_mex(i,4),pairwise_mex(i,5));

end
fclose(f);