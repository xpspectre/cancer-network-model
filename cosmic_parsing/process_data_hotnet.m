%% Process data for hotnet
%   Take TCGA breast cancer patients, compare to list of genes in hotnet
%   interaction database onlyl

load CosmicIncFus_samples_multi

load hotnetgenes

patient_data = {};
all_patient_genes = {}; % list of all patient mutations
acccepted_patient_genes = {}; % list of patient mutations in HPRD
for i = 1:size(samples_multi,1)
    name = samples_multi{i,1};
    tab = samples_multi{i,2};
    mutations = {};
    for j = 1:height(tab)
        gene = tab{j,1}{1};
        pos = find(strcmp(gene,hotnetgenes));
        if ~isempty(pos) % only take genes listed in hotnet
            mutations = [mutations; gene];
            acccepted_patient_genes = [acccepted_patient_genes; gene];
        end
        all_patient_genes = [all_patient_genes; gene];
    end
    patient_data = [patient_data; {name, mutations}];
end

save patient_data_multi patient_data

all_patient_genes = unique(all_patient_genes);
save all_patient_genes all_patient_genes

acccepted_patient_genes = unique(acccepted_patient_genes);
save acccepted_patient_genes all_patient_genes

%% Convert patient data matrix into format for hotnet
f = fopen('patient_mutations_multi.txt','w');
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