function MZ = pacejkaMZ(P, L, FZ, IA, alpha, FYexp)

%4.E31
dfz = (FZ - P(1)) ./ P(1);

%4.E32
Sht = (P(2) + P(3) .* dfz) .* L(1);
alphaT = alpha + Sht;

%4.E33
Ct = P(4) .* L(2);

%4.E34
Dt = (P(5) + P(6) .* dfz) .* (FZ ./ P(1)) .* L(3);

%4.E35
Et = (P(8)+ P(9) .* dfz) .* (1 - P(10) .* sign(alphaT)) .* L(5);

%4.E36-4.E38
Bt = (P(7) .* (1 + P(11) .* dfz)) ./ (Ct .* Dt) .* L(4);

%4.E39-4.E40
x_t = Bt .* alphaT;
t = Dt .* cos(Ct .* atan(x_t - Et .* (x_t - atan(x_t))));

%4.E41
Cr = P(12) .* L(7);

%4.E42
Dr = (P(13) + P(14) .* dfz) .* FZ .* L(6);

%4.E43-44
Br = P(15) ./ Cr .* L(8);

%4.E45-47
Mzr = Dr .* cos(Cr .* atan(Br .* alpha));

%Final
MZ = -t .* FYexp + Mzr;

end