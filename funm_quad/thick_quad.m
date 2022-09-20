function [ ell,U,T,D ] = thick_quad( U,T,D,ell )
%THICK_QUAD returns reordered Schur form of H
%   [U,T] = schur(H) and D = ordeig(H)
%   return reordered Schur form and eigenvalues such that
%   the ell "wanted" ones occur first
%
%   This file is part of the FUNM_QUAD code described in 
%
%   A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%   stable Arnoldi restarts for matrix functions based on quadrature,
%   SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%


[ignore,jj] = sort(abs(real(D)),1,'ascend');
ind = 0*D;
ind(jj(1:ell)) = 1;
[U,T] = ordschur(U,T,ind);
D = ordeig(T);

% do not break conjugate eigenvalues in the real Schur form
if(length(D)>ell && abs(conj(D(ell+1))-D(ell))<100*eps*abs(D(ell)) ),
    ell = ell + 1;
end;