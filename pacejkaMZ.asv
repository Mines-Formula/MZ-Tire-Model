function MZ0 = pacejkaMZ(P, L, FZ, IA, alpha, FY);

% inputs 24 Q coefficients,
% Estimated coefficients are as follows (found in Appen:
%Q = [0.007, -0.002, 0.147, 0.004, 8.964, -1.106, -0.842, 0, -0.227, 1.180, 0.1, -0.001, 0.007, 13.05, -1.609, -0.359, 0, 0.174, -0.896, 0, -0.008, 0, -0.296, -0.009];
%L = [1, 1, 1, 1, 1, 1, 1, 1];



t0 = Dt * cos(Ct * atan((Bt * alpha_t) - Et * ((Bt * alpha_t) - atan(Bt * alpha_t)))) * (Vcx / V_prime_c); % 4.E33

MZ0_prime = -t0 * FY; % 4.E32

MZ0 = MZ0_prime + Mzr0; % 4.E31