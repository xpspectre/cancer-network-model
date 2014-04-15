%% Load adjacency matrix
%   Should just contain 1's and 0's
load adjacency_interactions

num_proteins = size(A,1);

% Visualize interactions
% figure
% spy(A)
% xlabel('Protein')
% title('Adjacency Matrix')

%% Make degree matrix
d = sum(A);
D = spdiags(d(:),0,num_proteins,num_proteins);

%% Make Laplacian matrix
%   A is a multigraph (has loops but no multiedges) but doesn't affect calculation
L = D - A;

save master_network_matrics A A_index D L