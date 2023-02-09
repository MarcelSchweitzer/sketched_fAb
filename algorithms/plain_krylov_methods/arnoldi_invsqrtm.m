function [x, rel_err, err_est, iter] = arnoldi_invsqrtm(A, b, m, tol_stop, exact)

if nargin < 3
    m = 50;
end
if nargin < 4
    tol_stop = 1e-5;
end
if nargin < 5 || isempty(exact)
    exact = [];
    rel_err = [];
end

err_est = zeros(m,1);
iter = m;

N = size(b,1);
d = 20; % how often to check the error estimate

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
        coeffs = sqrtm(H(1:j,1:j))\(nb*e1);

        if j == d
            coeffs_old = coeffs;
        else
            err_est(j) = norm(coeffs-[coeffs_old; zeros(d,1)])/norm(coeffs);
            coeffs_old = coeffs;
            if err_est(j) < tol_stop
                iter = j;
                break
            end
        end
    end
end

if size(V,2) > length(coeffs)
    coeffs(size(V,2),1) = 0;
end
x = V*coeffs;
if ~isempty(exact)
    rel_err = norm(x-exact)/norm(exact);
end