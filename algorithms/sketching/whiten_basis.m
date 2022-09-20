function [SV, SAV, Rw] = whiten_basis(SV, SAV)

[Q,Rw] = qr(SV,0);
SV = Q;
SAV = SAV/Rw;