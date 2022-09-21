function hS = setup_sketching_handle(N,s)

E = spdiags(2*round(rand(N,1))-1,0,N,N); % Rademacher
D = speye(N); D = D(randperm(N,s),:);
hS = @(X) D*dct(E*X)/sqrt(s/N);
