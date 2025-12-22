function [Mz_pred, Kza0, Kzg0] = pacejka_MZ(P, Fz, alpha, gamma, Fy0, Kya0, R0)

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

df = Fz ./ max(Fz);
gamma_star = gamma;

SHt = qHz1 + qHz2 .* df.^2 + (qHz3 + qHz4 .* df.^2) .* gamma_star;

alpha_t alpha + SHt;

Dt = (Dt1 .* Fz + Dt2);
t0 = Dt .* cos(Ct.*atan(Bt .* alpha_t - Et .* (Bt .* alpha_t - atan(Bt .* alpha_t)))) .* cos(alpha);