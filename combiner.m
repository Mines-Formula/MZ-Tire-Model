clc, clearvars, clear all

dataFolder = 'data/Outputs/R20Round4&5';

filteredFiles = dir(fullfile(dataFolder, 'R20MZ_FZ_*_IA-4_filtered.csv'));

combinedTable = table();

for i = 1:numel(filteredFiles)
    curFile = fullfile(filteredFiles(i).folder, filteredFiles(i).name);
    disp("Adding " + filteredFiles(i).name);

    tempTable = readtable(curFile);

    [~, name, ~] = fileparts(filteredFiles(i).name);
    tempTable.SourceFile = repmat({name}, height(tempTable), 1);

    combinedTable = [combinedTable; tempTable];
end

combinedTable = sortrows(combinedTable, 'NormalForce');

outPath = fullfile(dataFolder, 'R20MZ_combined_IA-4_filtered.csv');
writetable(combinedTable, outPath);

disp("Combined filtered data saved as: " + outPath);
