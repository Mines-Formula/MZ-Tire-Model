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

%Normalized Load and camber
df = Fz ./ max(Fz); % normalized load
gamma_star = gamma; % normalized camber

%4.E35
SHt = qHz1 + qHz2 .* df.^2 + (qHz3 + qHz4 .* df.^2) .* gamma_star;

%4.E37
alpha_t alpha + SHt;

%4.E33
Dt = (Dt1 .* Fz + Dt2);
t0 = Dt .* cos(Ct.*atan(Bt .* alpha_t - Et .* (Bt .* alpha_t - atan(Bt .* alpha_t)))) .* cos(alpha);

%4.E32
Mz0_prime = -t0 .* Fy0;

%4.E36
alpha_r = alpha;
Dr = (Dr1 .* Fz + Dr2);
Mzr0 = Dr .* cos(Cr .* atan(Br .* alpha_r));

%4.E31
Mz_pred = Mz0_prime + Mzr0;

%4.E48
Kza0 = Dt .* Kya0;

%4.E49
Kzg0 = Fz .* R0 .* (qDz8 + qDz9 .* df.^2);

end