function [Mz_pred, Kza0, Kzg0] = pacejka_MZ(P, Fz, alpha, gamma, Fy0, Kya0, R0)

%Parameters
Bt = P(1);
Ct = P(2);
Dt1 = P(3);
Dt2 = P(4);
Et = P(5);
Br = P(6);
Cr = P(7);
Dr1 = P(8);
Dr2 = P(9);
qHz1 = P(10);
qHz2 = P(11);
qHz3 = P(12);
qHz4 = P(13);
qDz8 = P(14);
qDz9 = P(15);
qDz3 = P(16);
qDz4 = P(17);

%Normalized Load and camber
df = Fz ./ max(Fz); % normalized load
gamma_star = gamma; % normalized camber

%4.E35
SHt = qHz1 + qHz2 .* df.^2 + (qHz3 + qHz4 .* df.^2) .* gamma_star;

%4.E37
alpha_r = alpha_star + SHt; % apparently also equals alpha_f? Not sure what that is for

%4.E34
alpha_t = alpha_star + SHt

%4.E43
Dt0 = (Dt1 .* Fz + Dt2);
zeta5 = 1;
Dt = Dt0 .* (1 + qDz3 .* gamma_star + qDz4 .* gamma_star.^2) .* zeta5;

%4.E33
t0 = Dt .* cos(Ct .* atan(Bt.*alpha_t - Et.*(Bt.*alpha_t - atan(Bt.*alpha_t)))) .* cos(alpha);

%4.E32
Mz0_prime = -t0 .* Fy0;

%4.E36
Dr = (Dr1 .* Fz + Dr2); % This one needs to be double checked
Mzr0 = Dr .* cos(Cr .* atan(Br .* alpha_r)); % This is correct, just trying to figure out definitions for the parameters.

%4.E31
Mz0 = Mz0_prime + Mzr0;

%4.E48
Kza0 = Dt .* Kya0;

%4.E49
Kzg0 = Fz .* R0 .* (qDz8 + qDz9 .* df.^2);

%4.E39
K_prime_yAlpha = Kya + epsilon_k; % Not sure what epsilon+k is, it just shows up and then is never used again

%4.E38
SHf = SHy + SVy / K_prime_yAlpha;

%4.E40
Bt = (qBz1 + qBz2 .* fz + qBz3 + fz.^2) .* (1.+qBz4.*gamma_star + qBz5 .* abs(gamma_star)) .* gamma_Ky_alpha ./ gamma_mu_y; %Must be greater than zero

%4.E41
Ct = qCz1; % Must be greater than zero





end