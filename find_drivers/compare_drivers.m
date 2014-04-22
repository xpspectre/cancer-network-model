load all_mutations
load driver_mutations

matched_drivers = {};
for i = 1:length(all_mutations)
    pos = find(strcmp(all_mutations{i},drivers));
    if ~isempty(pos)
        matched_drivers = [matched_drivers; all_mutations{i}];
    end
end

f = fopen('matched_drivers.txt','w');
for i = 1:size(matched_drivers,1)
    fprintf(f,'%s\n',matched_drivers{i});
end
fclose(f);