function [ v1,H,h,breakdown,accuracy_flag ] = lanczos( A,m,H,s,param )
%LANCZOS extend a given Lanczos decomposition (V_big,H) of dimension s
%  to dimension m. This file has been adapted from the FUNM_KRYL code
%  described in
%
%   M. Afanasjew, M. Eiermann, O. G. Ernst, and S. G\"{u}ttel (2008):
%   Implementation of a Restarted Krylov Subspace Method for the Evaluation
%   of Matrix Functions, Linear Algebra Appl., 429:2293--2314.
%
%  It is now part of the FUNM_QUAD code described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%

accuracy_flag = 0;
fm = 0;
tol = param.tol;
if param.max_restarts == 1
    if strcmp(param.function,'invSqrt')
        fm = @(X) inv(sqrtm(X));
    end
    if strcmp(param.function,'log')
        fm = @(X) logm(X);
    end
    if strcmp(param.function,'exp')
        fm = @(X) expm(X);
    end
end

breakdown = 0;   % no checking
global V_big
H(m+1,m) = 0;

% orthog first mvp against ALL previous

v0 = V_big(:,s);
if isnumeric(A),
    w = A*v0;
else
    w = A(v0);
end;

for j = 1:s,
    ip = param.inner_product(w,V_big(:,j));
    H(j,s) = H(j,s) + ip(1);
    w = w - V_big(:,j)*ip(1);
end;

H(s+1,s) = sqrt(param.inner_product(w,w));
v1 = (1/H(s+1,s))*w;

% from here on: Lanczos with v0 and v1
for k = s+1:m,
    H(k-1,k) = H(k,k-1);
    V_big(:,k) = v1;
    if isnumeric(A),
        w = A*v1;
    else
        w = A(v1);
    end;
    w = w - H(k-1,k)*v0;
    H(k,k) = param.inner_product(w,v1);
    w = w - H(k,k)*v1;
    H(k+1,k) = sqrt(param.inner_product(w,w));
    if abs(H(k+1,k)) < k*eps*norm(H(1:k+1,k)),
        breakdown = k;
        break
    end
    v0 = v1;
    v1 = w*(1/H(k+1,k));
    if param.max_restarts == 1 && (~mod(k,10) && k >= 20),
        if isa(fm,'function_handle')
            c = fm(H(1:k,1:k))*eye(k,1);
        else
            [WW,DD] = eig(H(1:k,1:k));
            ee = unit(1,k);
            c = zeros(size(ee));
            for j = 1:k
                active_nodes = diag(DD);
                subdiag = diag(H(1:k+1,1:k),-1);
                fun = @(t) param.function(DD(j,j),t) .* evalnodal(t, active_nodes, subdiag);
                I = myintegral(fun,-inf,0,'AbsTol',tol,'RelTol',tol);
                c(j) = I;
            end
            c = (WW*spdiags(c,0,k,k)/WW)*ee;
        end
        if norm(c(k-9:k)) < norm(c)*param.stopping_accuracy/2,
            accuracy_flag = 1;
            breakdown = k;
            m = k;
            break
        end
    end
end;

h = H(m+1,m);
H = H(1:m,1:m);
