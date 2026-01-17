function MZ0 = pacejkaMZ(P, L, FZ, IA, alpha, FY);

% inputs 24 P coefficients,
% Estimated coefficients are as follows:
%P = [0.007, -0.002, 0.147, 0.004, 8.964, -1.106, -0.842, 0, -0.227, 1.180, 0.1, -0.001, 0.007, 13.05, -1.609, -0.359, 0, 0.174, -0.896, 0, -0.008, 0, -0.296, -0.009];
%L = [1, 1, 1, 1, 1, 1, 1, 1];



MZ0 = MZ0_prime + Mzr0; % 4.E31