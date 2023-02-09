function [err, diff_f] = sfom_closed_eval_error_and_diff(V,SV,SAV,A,Vfull,Sv,R,ex,f,num_it)

err = zeros(num_it,1);
diff_f = zeros(num_it,1);
for m = 1:num_it, m
    SVm = SV(:,1:m);
    M = SVm'*SVm;

    coeffs = M\(f( (SVm'*SAV(:,1:m))/M )*(SVm'*Sv));
    appr = V(:,1:m)*(R(1:m,1:m)\coeffs);
    err(m) = norm(appr - ex);
    f_sk =  R(1:m,1:m)\(M\(f((SVm'*SAV(:,1:m))/M)));
    f_tr = f(Vfull(:,1:m)'*A*Vfull(:,1:m));
    diff_f(m) = norm(f_sk-f_tr)/norm(f_tr);
end