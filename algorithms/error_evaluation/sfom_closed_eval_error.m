function err = sfom_closed_eval_error(V,SV,SAV,Sv,R,ex,f,num_it)

err = zeros(num_it,1);
for m = 1:num_it, m
    SVm = SV(:,1:m);
    M = SVm'*SVm;

    coeffs = M\(f( (SVm'*SAV(:,1:m))/M )*(SVm'*Sv));
    appr = V(:,1:m)*(R(1:m,1:m)\coeffs);
    err(m) = norm(appr - ex);
end