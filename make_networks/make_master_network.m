%% Load data
% Load list of mutated proteins, COSMIC/TCGA, breast cancer
load gene_id_list

% Load HPRD binary protein interactions
load interactions

%% Preprocess mutation data
%   Only keep proteins that show up in HPRD database (whether or not they interact with another member)
% Get list of all proteins in HPRD database
interacting_proteins = {};
for i = 1:height(interactions)
    p1 = interactions{i,1}{1};
    p2 = interactions{i,4}{1};
    if ~ismember(p1,interacting_proteins)
        interacting_proteins = [interacting_proteins; p1];
    end
    if ~ismember(p2,interacting_proteins)
        interacting_proteins = [interacting_proteins; p2];
    end
end
interacting_proteins = sort(interacting_proteins);
save interacting_proteins interacting_proteins

% Get list of proteins with mutations from TCGA breast cancer data
mutated_proteins = mutated_gene_id_pairs(:,2);
save mutated_proteins mutated_proteins

% Get subset of mutated proteins that appear in HPRD database
network_proteins = {};
for i = 1:length(mutated_proteins)
    if ismember(mutated_proteins{i},interacting_proteins)
        network_proteins = [network_proteins; mutated_proteins{i}];
    end
end
network_proteins = sort(network_proteins);
save network_proteins network_proteins

%% Make matrix of interactions A
%   Sparse matrix
%   1 for interaction, 0 otherwise
num_interacting_proteins = length(network_proteins);
idx_i = [];
idx_j = [];
for i = 1:height(interactions)
    
    % Get interacting pair
    p1 = interactions{i,1}{1};
    p2 = interactions{i,4}{1};
    
    % Only consider pairs where both proteins are network proteins
    p1pos = find(strcmp(p1,network_proteins));
    p2pos = find(strcmp(p2,network_proteins));
    if isempty(p1pos) || isempty(p2pos) % interaction not found
        continue
    elseif strcmp(p1,p2) % only add once to diagonal for self-interaction
        idx_i = [idx_i; p1pos];
        idx_j = [idx_j; p2pos];
    else % interacting pair found
        idx_i = [idx_i; p1pos; p2pos];
        idx_j = [idx_j; p2pos; p1pos];
    end
    
end

% Fill in sparse matrix
A = sparse(idx_i, idx_j, 1, num_interacting_proteins,num_interacting_proteins);

% Index of names of proteins for A
A_index = network_proteins;

num_interactions = length(idx_i)/2;

save adjacency_interactions A A_index


