function err = sgmres_invsqrt_quad_eval_error(V,SV,SAV,Sv,R,ex,num_it,c,z)

err = zeros(num_it,1);
ell = length(c);
for m = 1:num_it, m
    h = 0;
    SVm = SV(:,1:m);
    SAVm = SAV(:,1:m);
    for j = 1:ell
        X = (-(1-z(j))*SVm - (1+z(j))*SAVm);
        h = h + c(j)*( (X'*X)\(X'*Sv) );
    end
    h = -2/pi*h;
    appr = V(:,1:m)*(R(1:m,1:m)\h);
    err(m) = norm(appr - ex);
end