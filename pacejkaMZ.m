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

end