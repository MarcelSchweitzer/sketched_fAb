function err = best_approx_eval_error(V,ex,num_it)

err = zeros(num_it,1);
coeffs = V'*ex;
for m = 1:num_it
    appr = V(:,1:m)*coeffs(1:m);
    err(m) = norm(appr - ex);
end