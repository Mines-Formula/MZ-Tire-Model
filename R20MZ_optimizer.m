clc, clearvars, close all

data0 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-0_filtered.csv');
data2 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-2_filtered.csv');
data4 = readtable('data/Outputs/R20Round4&5/R20MZ_combined_IA-4_filtered.csv');
data  = [data0; data2; data4];
%data = data2;

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
    baseTag = sprintf('R20MZ_FZ_%d', FZBins(i));

    idx = contains(source, [baseTag '_IA-0_filtered']) | contains(source, [baseTag '_IA-2_filtered']) | contains(source, [baseTag '_IA-4_filtered']);

    MzAll = [MzAll; MzExperimental(idx)];
    FYAll = [FYAll; FYExperimental(idx)];
    FZAll = [FZAll; FZBins(i) * ones(sum(idx),1)];
    alphaAll = [alphaAll; alpha(idx)];
    IAAll = [IAAll; IA(idx)];
end

alphaAll = deg2rad(alphaAll);
IAAll = deg2rad(IAAll);

P_fixed = [250, 1.4, 2.4, -0.25, 3, -0.1, -1.5, 0, 0, -30.5, 1.15, 1, 0, 0, -0.128, 0, 0, 0, 1.43];
L = ones(1,8);

% Initial Q guess from Appendix 3
Q0 = [0.00276055, 0.00667048, 0.339315, -1.15563, 9.07658, -8.19445, -11.5486, -20.4342, 21.4985, 1.44589, 0.0184206, 0.0141927, 2.24868, -20.2342, 0.337115, -0.984545, -1.95786, 0.122852, -19.6941, -3.55172e-05, -0.000452427, 0.00165919, 0.295113, -0.375074];

objFun = @(Q) MzExperimental_all(P_fixed, Q, L, FZAll, IAAll, alphaAll, FYAll, MzAll);

options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 30000, 'TolFun', 1e-8, 'TolX', 1e-8);

lb = -Inf(size(Q0));
ub = Inf(size(Q0));

QOptimization = lsqnonlin(objFun, Q0, lb, ub, options);

disp('Optimized Q (Mz) Parameters:');
disp(QOptimization);

MzPredicted = pacejkaMZ(P_fixed, QOptimization, L, FZAll, IAAll, alphaAll, FYAll);

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

        [alphaSort, ord] = sort(alphaAll(idx));
        MzFit = pacejkaMZ(P_fixed, QOptimization, L, FZAll(idx), IAAll(idx), alphaSort, FYAll(idx));

        plot(rad2deg(alphaSort), MzFit, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Fit FZ=%d', FZBins(i)));
    end

    title(sprintf('Pacejka Mz Fit - IA = %d (RMSE = %.3f)', IAValues(k), rmse));
    xlabel('Slip Angle [deg]');
    ylabel('Aligning Torque Mz');
    legend('Location', 'best');
end

fprintf('[');
fprintf('%g, ', QOptimization(1:end-1));
fprintf('%g]\n', QOptimization(end));

function err = MzExperimental_all(P, Q, L, FZ, IA, alpha, FYexp, MzExp)
    MzModel = pacejkaMZ(P, Q, L, FZ, IA, alpha, FYexp);
    err = MzModel - MzExp;
end
