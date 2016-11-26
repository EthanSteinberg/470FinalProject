% Compute weights for every tree topology for the migration sample.
% The algorithm is to split up the SNP data into chunks and then use phylip
% to determine the topology of each chunk. Then you simply counts how many
% times each topology shows up and that is the weight.

data = tdfread('migrationSNP',' ');

trees = {};

for offset=1:100:size(data.DNA,2)
    fileID = fopen('infile','w');
    endIndex = min(size(data.DNA,2), offset+99);
    fprintf(fileID,'%d %d\n', size(data.SAMPLE, 1), endIndex - offset+1);
    for i=1:size(data.SAMPLE, 1)
        fprintf(fileID, '%s %s\n', data.SAMPLE(i,:), data.DNA(i, offset:endIndex));
    end
    fclose(fileID);

    !phylip dnaml < phylipInput > /dev/null
    
    outTreeID = fopen('outtree','r');
    trees{length(trees) + 1} = fscanf(outTreeID, '%s');
    fclose(outTreeID);
end

% Remove the lengths and RAT outgroup
for i=1:length(trees)
    trees{i} = regexprep(trees{i},':\d.\d\d\d\d\d','');
    trees{i} = regexprep(trees{i},',RAT','');
    trees{i} = trees{i}(1:end-1);
end

%% Count every tree
counts = containers.Map;
for i=1:length(trees)
    if ~isKey(counts, trees{i})
        counts(trees{i}) = 1;
    else
        counts(trees{i}) = counts(trees{i}) + 1;
    end
end


%% Output the weights

weightsID = fopen('weights','w');
fprintf(weightsID, 'TREE\tWEIGHT\n');

countsKeys = counts.keys;
countsValues = counts.values;

for i=1:length(countsKeys)
    treeString = countsKeys{i};
    weight = countsValues{i} / length(trees);
    
    fprintf(weightsID,'%s\t%f\n', treeString, weight);
end

fclose(weightsID);
