function MZ = pacejkaMZ(P, L, FZ, IA, alpha, FY)

dfz = (FZ - P(1)) ./ P(1); %4.E31
SHt = (P(2) + P(3).*dfz) .* L(1); %4.E32
alpha_t = alpha + SHt;
Ct = P(4) .* L(2); %4.E33
Dt = (P(5) + P(6).*dfz) .* (FZ ./ P(1)) .* L(3); %4.E34
Et = (P(8) + P(9).*dfz) .* (1 - (P(10) + P(11).*IA).*sign(alpha_t)) .* L(5); %4.E35
Bt = (P(7) .* (1 + P(12).*dfz)) ./ (Ct .* Dt) .* L(4); %4.E36-38
xt = Bt .* alpha_t %4.E39-4.E40
t = Dt .* cos( Ct .* atan(xt - Et .* (xt - atan(xt))));
Cr = P(13) .* L(6); %4.E41
Dr = (P(14) + P(15).*dfz) .* FZ .* L(7); %4.E42
Br = P(16) ./ Cr .* L(8); %4.E43-4.E44
Mzr = Dr .* cos(Cr .* atan(Br .* alpha)); %r.E45-47

MZ = -t .* FY + Mzr;

end