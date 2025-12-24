clc, clearvars, close all

data0 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-0_filtered.csv');
data2 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-0_filtered.csv');
data4 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-0_filtered.csv');
data  = [data0; data2; data4];

%These are different than the FY optimizer so pay attention
IA = data{:,3};
FZRaw = data{:,4};
alpha = data{:,5};
MzExperimental = data{:,7};
FYExperimental = data{:,8};
source = data{:,10};

FZBins = [50, 100, 150, 200, 250];

MzAll = [];
FYAll = [];
FZAll = [];
alphaAll = [];
IAAll = [];

for i = 1:numel(FZBins)
    baseTag = sprintf('R20_FZ_%d', FZBins(i));

    idx = contains(source, [baseTag '_IA-0_filtered']) | contains(source, [baseTag '_IA-2_filtered']) | contains(source, [baseTag '_IA-4_filtered']);

    MzAll = [MzAll; MzExperimental(idx)];
    FYAll = [FYAll; FYExperimental(idx)];
    FZAll = [FZAll; FZBins(i) * ones(sum(idx),1)];
    alphaAll = [alphaAll; alpha(idx)];
    IAAll = [IAAll; IA(idx)];
end

alphaAll = deg2rad(alphaAll);
IAAll = deg2rad(IAAll);

P0 = zeros(1, 15);
P0(1) = 250; % FZ0
P0(4) = 1.2; % Ct
P0(5) = 0.15; % Dt
P0(12) = 1.0; % Cr
P0(13) = 0.01; % Dr

L = ones(1,8);

objFun = @(P) MzExperimental_all(P, L, FZAll, IAAll, alphaAll, FYAll, MzAll);
options = optimoption('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 30000, 'TolFun', 1e-8, 'TolX', 1e-8);

lb = -Inf(size(P0));
ub = Inf(size(P0));

lb(1) = 250;
ub(1) = 250;

POptimization = lsqnonlin(objFun, P0, lb, ub, options);
disp('Optimized Mz Parameters:');
disp(POptimization);

MzPredicted = pacejkaMZ(POptimization, L, FZAll, IAAll, alphaAll, FYAll);

rmse = sqrt(mean((MzAll - MzPredicted).^2));
fprintf('Overall RMSE (Mz): %.4f Nm\n', rmse);

IAValues = [0, 2, 4];
colors = lines(numel(FZBins));

for k = 1:numel(IAValues)
    figure;
    hold on;
    grid on;

    thisIA = deg2rad(IAValues(k));
    idxIA = abs(IAAll - thisIA) < 1e-6;

    for i = 1:numel(FZBins)
        idx = idxIA & (FZAll == FZBins(i));

        scatter(rad2deg(alphaAll(idx)), MzAll(idx), 6, colors(i,:), 'filled', 'DisplayName', sprintf('Exp FZ=%d', FZBins(i)));

        [alphaSort, ord] = sort(alphaAll(id));
        MzFit = pacejkaMz(POptimization, L, FZAll(idx), IAAll(idx), alphaSort, FYAll(idx));

        plot(rad2deg(alphaSort), MZFit, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Fit FZ=%d', FZBins(i)));
    end

    title(sprintf('Pacejka Mz Fit - IA = %d (RMSE = %.3f)', IAValues(k), rmse));
    xlabel('Slip Angle [deg]');
    ylabel('Aligning Torque Mz');
    legend('Location', 'best');
end

fprintf('[');
fprintf('%g, ', POtimization(1:end-1));
fprint('%g]\n', POptimization(end));

