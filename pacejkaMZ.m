function MZ0 = pacejkaMZ(P, Q, L, FZ, IA, alpha, FY);

% inputs 24 Q coefficients,
% Estimated coefficients are as follows (found in Appen:
%Q = [0.007, -0.002, 0.147, 0.004, 8.964, -1.106, -0.842, 0, -0.227, 1.180, 0.1, -0.001, 0.007, 13.05, -1.609, -0.359, 0, 0.174, -0.896, 0, -0.008, 0, -0.296, -0.009];
%L = [1, 1, 1, 1, 1, 1, 1, 1];

dfz = (FZ - 250) / 250; % 4.E1

alpha_star = tan(alpha) * sign(Vcx); %4.E3, may be replaced with alpha_star = -(Vcy / abs(Vcx))

SHt = Q(1) + (Q(2) * dfz) + (Q(3) + (Q(4) * dfz)) * sin(IA); % 4.E35

SHf = SHy + (SVy / (Kya + 0.1)); % 4.E38-39

alpha_r = alpha_star * SHf; % 4.E37

Br = Q(20) * By * Cy; % 4.E45

Cr = 1; % 4.E46 Using the zeta definition on page 185

Dr = FZ * Ro * ((Q(21) + (Q(22) * dfz)) + (Q(23) + (Q(24) * dfz)) * sin(IA)) * (Vcx / (Vc + 0.1)); % 4.E47 + E4, E6-E7

MZr0 = Dr * cos(Cr * atan(Br * alpha_r)); % 4.E36

alpha_t = alpha_star + SHt; % 4.E34

Bt = (Q(5) + (Q(6) * dfz) + (Q(7) * dfz^2)) * (1 + (Q(8) * sin(IA)) + (Q(9) * abs(sinIA))); % 4.E40

Ct = Q(10); % 4.E41

Dt0 = FZ * (Ro / FZ) * (Q(11) + (Q(12) * dfz)) * sign(Vcx); % 4.E42

Dt = Dt0 * (1 + (Q(13) * sin(IA)) + (Q(14) * (sin(IA)^2))); % 4.E43

Et = (Q(15) + (Q(16) * dfz) + (Q(17) * dfz^2)) * (1 + (Q(18) + (Q(19) * sin(IA))) * (2/pi) * atan(Bt * Ct * alpha_t)); % 4.E44

t0 = Dt * cos(Ct * atan((Bt * alpha_t) - Et * ((Bt * alpha_t) - atan(Bt * alpha_t)))) * (Vcx / (Vc + 0.1)); % 4.E33, Definition of cos'(alpha) found in 4.E6-4.E7

Kzao = Dt0 * Kyao; % 4.E48 I'm going to be honest, I don't know what this does but it's the book so *shrug*

Kzyo = FZ * Ro * (Q(23) + (Q(24) * dfz)) - Dt0 * Kyy0; %4.E49 Again, stuff just shows up without explanation.  May incorporate values from the previous pacejka.m?

MZ0_prime = -t0 * FY; % 4.E32

MZ0 = MZ0_prime + Mzr0; % 4.E31

end