%% Process data for hotnet
%   Take TCGA breast cancer patients, compare to list of genes in hotnet
%   interaction database onlyl

load CosmicIncFus_samples_tcga

load hotnetgenes

patient_data = {};
for i = 1:size(samples_tcga,1)
    name = samples_tcga{i,1};
    tab = samples_tcga{i,2};
    mutations = {};
    for j = 1:height(tab)
        gene = tab{j,1}{1};
        pos = find(strcmp(gene,hotnetgenes));
        if ~isempty(pos) % only take genes listed in hotnet
            mutations = [mutations; gene];
        end
    end
    patient_data = [patient_data; {name, mutations}];
end

save patient_data patient_data

%% Convert patient data matrix into format for hotnet
f = fopen('patient_mutations.txt','w');
for i = 1:size(patient_data,1)
    name = patient_data{i,1};
    fprintf(f,'%s\t',name);
    genes = patient_data{i,2};
    for j = 1:length(genes)
        gene = genes{j};
        fprintf(f,'%s\t',gene);
    end
    fprintf(f,'\n');
end
fclose(f);