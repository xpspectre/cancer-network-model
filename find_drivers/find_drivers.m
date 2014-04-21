%% Get networks from hotnet

f = fopen('components.txt','r');

tline = fgets(f);
networks = {};
while ischar(tline)
    C = strsplit(tline,'\t');
    networks = [networks;{C}];
    tline = fgets(f);
end

fclose(f);

%% Get mutations for patients
load patient_mutations

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
chis = [];
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
        joint_count = sum(net_mat(p1,:).*net_mat(p2,:)); % number of times mutation appears in both patients at once
        joint_prob = joint_count/num_patients;
        chi = joint_prob/(gene_probs(p1)*gene_probs(p2));
        chis = [chis; i, p1, p2, chi];
    end
end

