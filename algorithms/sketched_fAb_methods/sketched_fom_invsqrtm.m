function [x, rel_err, num_quad, err_est, iter] = sketched_fom_invsqrtm(A, b, m, s, trunc, tol_stop, tol_quad, use_quad, exact)

if nargin < 3
    m = 50;
end
if nargin < 4
    s = 2*m+1;
end
if nargin < 5
    trunc = 2;
end
if nargin < 6
    tol_stop = 1e-5;
end
if nargin < 7
    tol_quad = 1e-1*tol_stop;
end
if nargin < 8
    use_quad = false;
end
if nargin < 9 || isempty(exact)
    exact = [];
    rel_err = [];
end

num_quad = zeros(m,1);
err_est = zeros(m,1);
iter = m;

N = size(b,1);
d = 20; % how often to check the error estimate

D = speye(N);
D = D(randperm(N,s),:);
e = rand(N,1); e = (e >= 0.5) - (e < 0.5);
E = spdiags(e,0,N,N);
S = @(M) sqrt(N/s)*(D*(dct(E*M)));

SV = zeros(s,m);
SAV = zeros(s,m);

Sb = S(b);

V = zeros(N,m+1);
V(:,1) = b/norm(b);

H = zeros(m+1,m);

N1 = 32; N2 = 45;
current_eps = 0;

for j = 1:m
    if isnumeric(A)
        w = A*V(:,j);
    else
        w = A(V(:,j));
    end

    SV(:,j) = S(V(:,j));
    norm_SV2 = norm(SV(:,j))^2;
    ratio = 1/norm_SV2;
    if ratio < 1
        new_eps = 1-ratio;
    else
        new_eps = ratio-1;
    end
    if new_eps > current_eps
        current_eps = new_eps;
    end

    SAV(:,j) = S(w);

    H(1:j,j) = 0;
    for i = max([1,j-trunc+1]):j
        h = V(:,i)'*w;
        w = w - V(:,i)*h;
        H(i,j) = H(i,j) + h;
    end

    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j);

    if mod(j,d) == 0
        [Q,Rw] = qr(SV(:,1:j),0);
        SVw = Q;
        SAVw = SAV(:,1:j)/Rw;

        if use_quad
            quad_err = inf;
            first_time = true;
            while quad_err > tol_quad && N2 < 1000
                if first_time == true
                    weights1 = pi/N1*ones(1,N1);
                    t1 = zeros(1,N1); t2 = zeros(1,N2);
                    for ii = 1:N1
                        t1(ii) = cos((2*ii-1)/(2*N1) * pi);
                    end
                    f1 = zeros(j,1);
                    for ii = 1:N1
                        X = -(1-t1(ii))*SVw - (1+t1(ii))*SAVw;
                        f1 = f1 + weights1(ii)*( (X'*X)\(X'*Sb) );
                    end
                    f1 = -2/pi*f1;
                else
                    f1 = f2;
                end
                weights2 = pi/N2*ones(1,N2);
                for ii = 1:N2
                    t2(ii) = cos((2*ii-1)/(2*N2) * pi);
                end
                f2 = zeros(j,1);
                for ii = 1:N2
                      X = -(1-t2(ii))*SVw - (1+t2(ii))*SAVw;
                      f2 = f2 + weights2(ii)*((X'*X)\(X'*Sb));
                end
                f2 = -2/pi*f2;
                quad_err = norm(f1-f2)/norm(f2);
                num_quad(j) = N2;
                if quad_err > tol_quad
                    N1 = N2;
                    N2 = floor(sqrt(2)*N2);
                else
                    if first_time == true
                        N2 = N1;
                        N1 = ceil(N1/sqrt(2));
                    end
                end
                first_time = false;
            end
            if N2 > 1000
                warning(string(['Quadrature accuracy not reached. Attained accuracy is ', num2str(quad_err)]))
                N2 = N1;
                N1 = ceil(N1/sqrt(2));
            end
            coeffs = f2;
        else
            M = SVw'*SVw;
            coeffs = M\(sqrtm( (SVw'*SAVw)/M )\(SVw'*Sb));
        end

        if j == d
            coeffs_old = coeffs;
        else
            err_est(j) = 1/sqrt(1-current_eps)*norm(SV(:,1:j)*(coeffs-[coeffs_old; zeros(d,1)]))/norm(SV(:,1:j)*coeffs);
            coeffs_old = coeffs;
            if err_est(j) < tol_stop
                iter = j;
                break
            end
        end

    end

end

y = Rw\coeffs;
if size(V,2) > length(y)
    y(size(V,2),1) = 0;
end
x = V*y;
if ~isempty(exact)
    rel_err = norm(x-exact)/norm(exact);
end

