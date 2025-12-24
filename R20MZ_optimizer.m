clc, clearvars, close all

data0 = readtable('data/R20MZ_combined_filtered.csv');
data2 = readtable('data/R20MZ_IA-2_combined_filtered.csv');
data4 = readtable('data/R20MZ_IA-4_combined_filtered.csv');
data  = [data0; data2; data4];