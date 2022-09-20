function [ w,H,h,breakdown,accuracy_flag ] = arnoldi( A,m,H,s,param )
%ARNOLDI   Extend a given Arnoldi decomposition (V_big,H) of dimension s
%  to dimension m. This file has been adapted from the FUNM_KRYL code
%  described in 
%
%  M. Afanasjew, M. Eiermann, O. G. Ernst, and S. G\"{u}ttel (2008):
%  Implementation of a restarted Krylov subspace method for the evaluation
%  of matrix functions, Linear Algebra Appl., 429:2293--2314.
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

global V_big
H(m+1,m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;

for k = s:m,
    
    w = V_big(:,k);
    if isnumeric(A),
        w = A*w;
    else
        w = A(w);
    end
    
    sj = max([1,k-trunc+1]):k;
    if(k==s), sj = 1; end
    for r = 0:reo,
        for j = sj:k,
            ip = param.inner_product(w,V_big(:,j));
            H(j,k) = H(j,k) + ip(1);
            w = w - V_big(:,j)*ip(1);
        end
    end
    
    H(k+1,k) = sqrt(param.inner_product(w,w));
    
    if abs(H(k+1,k)) < k*eps*norm(H(1:k+1,k))
        breakdown = k;
        break
    end
    
    w = (1/H(k+1,k))*w;
    if k < m
        V_big(:,k+1) = w;
    end
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
    
end

h = H(m+1,m);
H = H(1:m,1:m);
