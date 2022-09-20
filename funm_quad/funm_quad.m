function [f,out,param] = funm_quad(A,b,param)
% FUNM_QUAD   Approximate f(A)*b by quadrature-based restarted Arnoldi
%
% This method is described in detail in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% and analyzed in 
%
% Parts of this code are based on the FUNM_KRYL method described in
%
%  M. Afanasjew, M. Eiermann, O. G. Ernst, and S. G\"{u}ttel (2008):
%  Implementation of a restarted Krylov subspace method for the evaluation
%  of matrix functions, Linear Algebra Appl., 429:2293--2314.
%
% Please refer to the demo files to see how this code is used.

if nargin < 3,
    param = struct;
end

[param,modified] = param_init_quad(param);

if modified && param.verbose,
    disp('Some of the parameters were inconsistent. FUNM_QUAD will now run with a corrected setting.');
end

fun_switch = 0;
if strcmp(param.function,'invSqrt'),
    fun_switch = 1;
end
if strcmp(param.function,'log')
    fun_switch = 2;
end
if strcmp(param.function,'exp')
    fun_switch = 3;
end
if isa(param.function,'function_handle')
    fun_switch = 4;
end

if ~fun_switch,
    error('FUNM_QUAD: param.function is neither a built-in function nor a function handle.');
end


if param.waitbar,
    hand = waitbar(0,'Please wait...');
end

beta = sqrt(param.inner_product(b,b));
v = (1/beta)*b;

n = length(b);
if param.restart_length >= n,
    if param.verbose >= 1,
        disp('Warning: Restart length larger than matrix dimension. Running Arnoldi/Lanczos without restarts.');
    end
    param.restart_length = n;
    param.max_restarts = 1;
end

m = param.restart_length;
ell = 0;    % thick restart param
subdiag = []; % subdiagonal entries of hessenberg matrices (for
% computing the norm of v_{m+1})

f = zeros(size(b));
if (param.H_full)
    H_full = [];
end

if (param.V_full)
    out.V_full = [];
end

out.stop_condition = '0 - maximal # of restarts exceeded';

%%

% allocate memory [ >= m+max(ell) ]
global V_big
alloc = param.restart_length + 20;
V_big = zeros(length(b),alloc);


N = 32; % initial number of quadrature points
if strcmp(param.function,'invSqrt')
    beta_transform = param.transformation_parameter;
end

% restart loop starts here
for k = 1:param.max_restarts,
    % check whether a stop condition is satisfied
    if str2double(out.stop_condition(1)),
        break
    end
    
    if param.verbose >= 2,
        disp(['---- Restart cycle ',num2str(k),' ----']);
    end
    
    % waitbar and cputime
    if param.waitbar,
        waitbar(k/param.max_restarts,hand);
    end
    out.time(k) = cputime;
    
    % compute A-invariant subspace of prev. cycle (thick restart)
    if (~isempty(param.thick) && k > 1),
        ell_prev = ell;
        [ ell,U,T,D ] = param.thick( U,T,D, param.number_thick );
        out.thick_replaced{k-1} = D(1:ell);
        out.thick_interpol{k-1} = D(ell+1:end);
        if ell,
            U = U(:,1:ell);
            V_big(:,1:ell) = V_big(:,1:m+ell_prev)*U;
            H_hat = U'*H*U;
            H = [H_hat; eta*U(end,:)];
        else
            H = [];
        end
    else
        H = [];
    end
    
    V_big(:,ell+1) = v;
    
    % compute/extend Krylov decomposition
    if param.hermitian,
        [ v,H,eta,breakdown, accuracy_flag ] = lanczos( A,m+ell,H,ell+1,param );
    else
        [ v,H,eta,breakdown, accuracy_flag ] = arnoldi( A,m+ell,H,ell+1,param );
    end
    
    if breakdown
        m = breakdown;
        if ~accuracy_flag
            if param.verbose >= 2,
                if param.hermitian == 0,
                    disp(['Arnoldi breakdown in cycle ', num2str(k),', iteration ',num2str(breakdown)]);
                else
                    disp(['Lanczos breakdown in cycle ', num2str(k),', iteration ',num2str(breakdown)]);
                end
            end
            H = H(1:breakdown+ell,1:breakdown+ell);
        end
    end
    
    out.thick_interpol{k} = eig(H);
    
    % store full Krylov basis?
    if param.V_full,
        out.V_full = [ out.V_full , V_big(:,1:m+ell) ];
    end
    
    % Schur form of H (for thick restarts)
    if (~isempty(param.thick)),
        if isreal(H),
            [U,T] = schur(H,'real');
        else
            [U,T] = schur(H);
        end
        D = ordeig(T);
    end
    
    % interpolation nodes currently active in f:
    %    k = 1 -> []
    %    k = 2 -> [ out1.thick_interpol{1} ; out1.thick_replaced{1} ]
    %    k = 3 -> [ out1.thick_interpol{1} ; out1.thick_interpol{2} ; out1.thick_replaced{2} ]
    %    ...
    % (note that the "replaced nodes" are going to be replaced in the next
    % sweep, respectively)
    
    active_nodes = [];
    for kk = 1:k-1,
        active_nodes = [ active_nodes ; out.thick_interpol{kk} ];
    end
    if ~isempty(param.thick) && k > 1,
        active_nodes = [ active_nodes ; out.thick_replaced{k-1} ];
    end
    
    if (param.H_full),
        H_full = blkdiag(H_full,H);
        if k>1,
            H_full(end-m+1,end-m-ell) = s;
        end
    end
    
    % approximate the restart function:
    % "h = e_k(H)*unit(ell+1,m+ell)"
    if k == 1 && fun_switch ~= 4 % in the first iteration g_1 = f and we use the "closed" form
        switch fun_switch
            case 1
                h2 = sqrtm(H)\unit(ell+1,m+ell);
            case 2
                h2 = logm(eye(m+ell)+H)*(H\unit(ell+1,m+ell));
            case 3
                h2 = expm(H)*unit(ell+1,m+ell);
        end
    else
        converged = 0;
        tol = param.tol;
        h1 = [];
        while ~converged,
            % compute eigenvalue decomposition of Hessenberg matrix
            [WW,DD] = eig(H);
            
            % approximate the restart function by Gauss-Jacobi quadrature
            if isempty(h1) && fun_switch ~= 4,
                N2 = N;
                if N > 2
                    N=floor(N/sqrt(2));
                end
                if mod(N,2) == 1
                    N = N-1;
                end
                
                switch fun_switch
                    case 1
                        % if f(z) = 1/sqrt(z), use Gauss-Jacobi quadrature
                        weights = pi/N*ones(1,N);
                        t = zeros(1,N);
                        for ii = 1:N
                            t(ii) = cos((2*ii-1)/(2*N) * pi);
                        end
                    case 2
                        %if f(z) = log(1+z)/z use Gauss-Legendre quadrature
                        b_gauss = .5./sqrt(1-(2*(1:N-1)).^(-2));
                        T_gauss = diag(b_gauss,1) + diag(b_gauss,-1);
                        % compute nodes and weights by solving a
                        % tridiagonal eigenproblem.
                        % WARNING: slow for large N
                        [V_gauss,D_gauss] = eig(T_gauss);
                        t = diag(D_gauss); [t,ind] = sort(t);
                        weights = 2*V_gauss(1,ind).^2;
                    case 3
                        % if f(z) = exp(z), use midpoint rule on parabolic
                        % Hankel contour
                        aa = max(1,max(real(active_nodes))+1);
                        bb = 1;
                        thetas = imag(active_nodes);
                        ccs = abs((active_nodes - aa - 1i*thetas)./thetas.^2);
                        cc = min(ccs)/5; % safety parameter
                        cc = min(cc,0.25);
                        phi = @(theta) aa + 1i*bb*theta - cc*theta.^2;
                        thetac = sqrt((aa-log(tol))/cc); % critical theta
                        theta = linspace(-thetac,thetac,N);
                        hh = theta(2)-theta(1);
                        z = phi(theta); % quad nodes
                        c = -hh/(2i*pi)*exp(z).*(1i*bb - 2*cc*theta);
                end
                
                % Evaluate the reciprocal of the nodal polynomial at the
                % quadrature points
                switch fun_switch
                    case 1
                        tt = -beta_transform*(1-t)./(1+t);
                    case 2
                        tt = -2./(t+1);
                    case 3
                        tt = z(1:N/2);
                end
                
                if isempty(param.thick)
                    rho_vec = evalnodal(tt, active_nodes, subdiag).';
                else
                    rho_vec = evalnodal(tt, active_nodes(1:end-length(out.thick_replaced{k-1})), subdiag(1:end-length(out.thick_replaced{k-1}))).';
                    rho_vec_replaced = evalnodal(tt, out.thick_replaced{k-1}, subdiag(end-length(out.thick_replaced{k-1})+1:end)).';
                    rho_vec = rho_vec .* rho_vec_replaced;
                end
                
                if param.hermitian % for Hermitian matrices, use diagonalization and scalar quadrature
                    ee = unit(ell+1,ell+m);
                    ww = WW\ee;
                    h1 = zeros(size(ee));
                    switch fun_switch
                        case 1
                            for j = 1:length(t)
                                h1 = h1 + weights(j)*rho_vec(j)*((-beta_transform*(1-t(j))*eye(ell+m)-DD*(1+t(j)))\ww);
                            end
                            h1 = -2*sqrt(beta_transform)/pi*h1;
                        case 2
                            for j = 1:length(t)
                                h1 = h1 + weights(j)*rho_vec(j)*((DD*(1+t(j))+2*eye(ell+m))\ww);
                            end
                        case 3
                            c = c(1:N/2).*rho_vec.';
                            h1 = zeros(m+ell,1);
                            for j = 1:N/2
                                h1 = h1 - c(j)*((z(j)*speye(size(DD))-DD)\ww);
                            end
                            h1 = 2*real(h1);
                    end
                    h1 = WW*h1;
                else % for non-Hermitian matrices, use matrix quadrature to avoid diagonalization
                    ee = unit(ell+1,ell+m);
                    h1 = zeros(size(ee));
                    switch fun_switch
                        case 1
                            for j = 1:length(t)
                                h1 = h1 + weights(j)*rho_vec(j)*((-beta_transform*(1-t(j))*eye(ell+m)-H*(1+t(j)))\ee);
                            end
                            h1 = -2*sqrt(beta_transform)/pi*h1;
                        case 2
                            for j = 1:length(t)
                                h1 = h1 + weights(j)*rho_vec(j)*((H*(1+t(j))+2*eye(m+ell))\ee);
                            end
                        case 3
                            c = c(1:N/2).*rho_vec.';
                            for j = 1:N/2
                                h1 = h1 - c(j)*((z(j)*speye(size(H))-H)\ee);
                            end
                            h1 = 2*real(h1);
                    end
                end
            end
            
            N2 = ceil(sqrt(2)*N);
            if mod(N2,2) == 1
                N2 = N2+1;
            end
            
            switch fun_switch
                case 1
                    % if f(z) = 1/sqrt(z), use Gauss-Jacobi quadrature
                    weights2 = pi/N2*ones(1,N2);
                    t2 = zeros(1,N2);
                    for ii = 1:N2
                        t2(ii) = cos((2*ii-1)/(2*N2) * pi);
                    end
                case 2
                    %if f(z) = log(1+z)/z use Gauss-Legendre quadrature
                    b_gauss2 = .5./sqrt(1-(2*(1:N2-1)).^(-2));
                    T_gauss2 = diag(b_gauss2,1) + diag(b_gauss2,-1);
                    % compute nodes and weights by solving a
                    % tridiagonal eigenproblem.
                    % WARNING: slow for large N
                    [V_gauss2,D_gauss2] = eig(T_gauss2);
                    t2 = diag(D_gauss2); [t2,ind] = sort(t2);
                    weights2 = 2*V_gauss2(1,ind).^2;
                case 3
                    aa = max(1,max(real(active_nodes))+1);
                    bb = 1;
                    thetas = imag(active_nodes);
                    ccs = abs((active_nodes - aa - 1i*thetas)./thetas.^2);
                    cc = min(ccs)/5; % safety parameter
                    cc = min(cc,0.25);
                    phi = @(theta) aa + 1i*bb*theta - cc*theta.^2;
                    thetac = sqrt((aa-log(tol))/cc); % critical theta
                    theta2 = linspace(-thetac,thetac,N2);
                    hh2 = theta2(2)-theta2(1);
                    z2 = phi(theta2); % quad nodes
                    c2 = -hh2/(2i*pi)*exp(z2).*(1i*bb - 2*cc*theta2);
            end
            
            % Evaluate the reciprocal of the nodal polynomial at the
            % quadrature points
            switch fun_switch
                case 1
                    tt = -beta_transform*(1-t2)./(1+t2);
                case 2
                    tt = -2./(t2+1);
                case 3
                    tt = z2(1:N2/2);
            end
            
            if fun_switch ~= 4
                if isempty(param.thick)
                    rho_vec2 = evalnodal(tt, active_nodes, subdiag).';
                else
                    rho_vec2 = evalnodal(tt, active_nodes(1:end-length(out.thick_replaced{k-1})), subdiag(1:end-length(out.thick_replaced{k-1}))).';
                    rho_vec_replaced2 = evalnodal(tt, out.thick_replaced{k-1}, subdiag(end-length(out.thick_replaced{k-1})+1:end)).';
                    rho_vec2 = rho_vec2 .* rho_vec_replaced2;
                end
            end
            
            if param.hermitian || fun_switch == 4   %for Hermitian matrices, use diagonalization and scalar quadrature
                ee = unit(ell+1,ell+m);
                ww = WW\ee;
                h2 = zeros(size(ee));
                switch fun_switch
                    case 1
                        for j = 1:length(t2)
                            h2 = h2 + weights2(j)*rho_vec2(j)*((-beta_transform*(1-t2(j))*eye(ell+m)-DD*(1+t2(j)))\ww);
                        end
                        h2 = -2*sqrt(beta_transform)/pi*h2;
                    case 2
                        for j = 1:length(t2)
                            h2 = h2 + weights2(j)*rho_vec2(j)*((DD*(1+t2(j))+2*eye(ell+m))\ww);
                        end
                    case 3
                        c2 = c2(1:N2/2).*rho_vec2.';
                        for j = 1:N2/2
                            h2 = h2 - c2(j)*((z2(j)*speye(size(DD))-DD)\ww);
                        end
                        h2 = 2*real(h2);
                    case 4
                        out.num_quadpoints(k) = 0;
                        max_err = 0;
                        for j = 1:m+ell
                            fun = @(t) param.function(DD(j,j),t) .* evalnodal(t, active_nodes, subdiag);
                            [I,errbnd,npts] = myintegral(fun,-inf,0,'AbsTol',tol,'RelTol',tol);
                            out.num_quadpoints(k) = max([npts out.num_quadpoints(k)]);
                            max_err = max(max_err,errbnd);
                            h2(j) = I;
                        end
                        if param.verbose >= 2,
                            disp([num2str(out.num_quadpoints(k)),' quadrature points were used. Error bound: ', num2str(max_err)])
                        end
                        h2 = (WW*spdiags(h2,0,m+ell,m+ell)/WW)*ee;
                        converged = 1;
                end
                if fun_switch ~= 4
                    h2 = WW*h2;
                end
            else %for non-Hermitian matrices, use matrix quadrature to avoid diagonalization
                ee = unit(ell+1,ell+m);
                h2 = zeros(size(ee));
                switch fun_switch
                    case 1
                        for j = 1:length(t2)
                            h2 = h2 + weights2(j)*rho_vec2(j)*((-beta_transform*(1-t2(j))*eye(ell+m)-H*(1+t2(j)))\ee);
                        end
                        h2 = -2*sqrt(beta_transform)/pi*h2;
                    case 2
                        for j = 1:length(t2)
                            h2 = h2 + weights2(j)*rho_vec2(j)*((H*(1+t2(j))+2*eye(m+ell))\ee);
                        end
                    case 3
                        c2 = c2(1:N2/2).*rho_vec2.';
                        for j = 1:N2/2
                            h2 = h2 - c2(j)*((z2(j)*speye(size(H))-H)\ee);
                        end
                        h2 = 2*real(h2);
                end
            end
            
            % Check if quadrature rule has converged
            if fun_switch ~= 4
                if norm(beta*(h2-h1))/norm(f) < tol
                    if param.verbose >= 2,
                        disp([num2str(N),' quadrature points were enough. Norm: ', num2str(norm(h2-h1)/norm(f))])
                    end
                    out.num_quadpoints(k) = N;
                    converged = 1;
                else
                    if param.verbose >= 2,
                        disp([num2str(N),' quadrature points were not enough. Trying ',num2str(N2),'. Norm: ', num2str(norm(h2-h1)/norm(f))])
                    end
                    h1 = h2;
                    N = N2;
                end
            end
        end
    end
    % end of restart function approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % workaround due to matlab 'bug' (copies large submatrices)
    
    h_big = beta*h2(1:m+ell,1);
    if size(V_big,2) > length(h_big),
        h_big(size(V_big,2),1) = 0;
    end
    % update Krylov approximation
    f = V_big*h_big + f;
    
    out.appr(:,k) = f;
    out.update(k) = beta*norm(h_big);  % norm of update
    
    
    % keep track of subdiagonal entries of Hessenberg matrix
    if m ~= 1
        subdiag = [ subdiag ; diag(H(end-m+1:end,end-m+1:end),-1) ; eta ];
    else
        subdiag = [ subdiag ; eta ];
    end
    s = eta;
    
    
    % get cpu-time
    out.time(k) = cputime - out.time(k);
    
    % check stopping condition
    if  stopcondition(out.update/norm(f) < param.stopping_accuracy) || accuracy_flag % stop by norm of update?
        out.stop_condition = '5 - norm of updates decayed below stopping accuracy';
    end
    
    if ~isempty(param.exact),
        % stop by absolute error?
        out.err(k) = norm(f - param.exact);
        
        if stopcondition(out.err(2:end)./out.err(1:end-1) >  param.min_decay),
            out.stop_condition = '3 - linear convergence rate of absolute error > min_decay';
        end
        if out.err(k) < param.stopping_accuracy,
            out.stop_condition = '1 - absolute error below stopping accuracy';
        end
    end
    
    if breakdown && ~accuracy_flag,
        out.stop_condition = '6 - Arnoldi/Lanczos breakdown';
    end
    
end
% restart loop ends here

if param.H_full,
    out.H_full = H_full;
end


if param.waitbar,
    close(hand);
end

clear V_big
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = stopcondition(v)
% stop condition
% checks whether the entries of v are of the form 1,...,1,0,...0,1,1

r = 0;

if length(v) < 2,
    return
end

if all(v(end-1:end)) && ~all(v),
    r = 1;
    return
end

end
