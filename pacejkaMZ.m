function MZ = pacejkaMZ(P, L, FZ, IA, alpha, FYexp)

dfz = (FZ - P(1)) ./ P(1);
Sht = (P(2) + P(3) .* dfz) .* L(1);
alphaT = alpha + Sht;
Ct = P(4) .* L(2);
Dt = (P(5) + P(6) .* dfz) .* (FZ ./ P(1)) .* L(3);
Et = (P(8)+ P(9) .* dfz) .* (1 - P(10) .* sign(alphaT)) .* L(5);
Bt = (P(7) .* (1 + P(11) .* dfz)) ./ (Ct .* Dt) .* L(4);
x_t = Bt .* alphaT;
t = Dt .* cos(Ct .* atan(x_t - Et .* (x_t - atan(x_t))));
Cr = P(12) .* L(7);
Dr = (P(13) + P(14) .* dfz) .* FZ .* L(6);
Br = P(15) ./ Cr .* L(8);
Mzr = Dr .* cos(Cr .* atan(Br .* alpha));
MZ = -t .* FYexp + Mzr;

end