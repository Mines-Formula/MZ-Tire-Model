clc, clearvars, close all

data0 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_filtered.csv');
data2 = readtable('data/Outputs/R20Round4&5/R20MZ_IA-2_combined_filtered.csv');
data4 = readtable('data/Outputs/R20Round4&5/R20MZ_IA-4_combined_filtered.csv');
data  = [data0; data2; data4];

%These are different than the FY optimizer so pay attention
IA = data{:,3};
FZRaw = data{:,4};
alpha = data{:,5};
MZ = data{:,7};
FYExperimental = data{:,8};
source = data{:,10};

FZBins = [50, 100, 150, 200, 250];

MZAll = [];
FYAll = [];
FZAll = [];
alphaAll = [];
IAAll = [];

for i = 1:numel(FZBins)
    baseTag = sprintf('R20_FZ_%d', FZBins(i));

    idx = contains(source, [baseTag '_filtered']) | contains(source, [baseTag '_IA-2_filtered']) | contains(source, [baseTag '_IA-4_filtered']);