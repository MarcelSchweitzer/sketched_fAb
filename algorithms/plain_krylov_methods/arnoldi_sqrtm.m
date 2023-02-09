function [x, err, iter] = arnoldi_sqrtm(A, b, m, tol_stop, exact)

if nargin < 3
    m = 50;
end
if nargin < 4
    tol_stop = 1e-5;
end
if nargin < 5 || isempty(exact)
    exact = [];
end

d = 20;
err = zeros(floor(m/d),1);
iter = m;


N = size(b,1);

V = zeros(N,m+1);
nb = norm(b);
V(:,1) = b/nb;

H = zeros(m+1,m);

for j = 1:m
    if isnumeric(A)
        w = A*V(:,j);
    else
        w = A(V(:,j));
    end

    H(1:j,j) = 0;
    for i = 1:j
        h = V(:,i)'*w;
        w = w - V(:,i)*h;
        H(i,j) = H(i,j) + h;
    end

    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j);

    if mod(j,d) == 0
        e1 = eye(j,1);
        coeffs = sqrtm(H(1:j,1:j))*(nb*e1);
        if size(V,2) > length(coeffs)
            coeffs(size(V,2),1) = 0;
        end
        x = V*coeffs;
    
        err(j/d) = norm(x-exact);
        if err(j/d) < tol_stop
            iter = j;
            break
        end
    end
end