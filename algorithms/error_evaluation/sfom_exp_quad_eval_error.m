function err = sfom_exp_quad_eval_error(V,SV,SAV,Sv,R,ex,num_it,c,z)

err = zeros(num_it,1);
ell = length(c);
for m = 1:num_it, m
    SVm = SV(:,1:m);
    SAVm = SAV(:,1:m);
    h = 0;
    for j = 1:ell/2
        X = SVm'*(z(j)*SVm - SAVm);
        h = h + c(j)*( (X'*X)\(X'*(SVm'*Sv)) ) ;
    end
    h = 2*real(h);
    appr = V(:,1:m)*(R(1:m,1:m)\h);
    err(m) = norm(appr - ex);
end