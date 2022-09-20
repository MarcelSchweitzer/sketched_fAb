%% Experiment from "Section 5.2 - Network example"
%
% This script is not optimized for best performance, as its primary
% purpose is generating plots that illustrate the convergence behavior
% of methods.

% Add paths etc.
addpath(genpath('algorithms'),'auxiliary','funm_quad');
mydefaults

% load matrix and solution
load('data/wiki-Vote.mat'); load('data/wiki-Vote-comp.mat'); % expm(-A)*v where v=ones(8797,1)

% initializations
A = -Problem.A; ee = -ee; N = size(A,1); 
v = ones(N,1);


% Number of Arnoldi iterations to perform
num_it = 50;

%% full Arnoldi and best approx
[~,~,~,Vfull] = bta(A,speye(N),v,inf(1,num_it),num_it,@(X) X);
err_best = best_approx_eval_error(Vfull,ex_expm,num_it);

%% truncated Arnoldi
for k = 2:4 % how many blocks to orthogonalise against
    figure(k-1)
    h0 = semilogy(err_best,'k');
    hold on

    hS = setup_sketching_handle(N,2*num_it); % s = 2*m_max
    [SV,SAV,~,Vtrunc] = bta(A,speye(N),v,inf(1,num_it),k,hS); % number of mat-vec products = num_it

    % whitening the basis
    [SV, SAV, Rw] = whiten_basis(SV, SAV);


    %% set up quadrature rule
    %  one to rule them all (determined for m_max)

    [Q,R] =  qr(SV,0);
    [X,Lam] = eig(Q'*SAV,R); Lam = diag(Lam);
    active_nodes = Lam;

    tol = 1e-16;
    ell = 100;
    [z,c] = exp_quad(active_nodes,tol,ell);

    if k == 4
        % Plot Ritz values, eigenvalues and integration contour
        figure(4)
        plot(real(ee),imag(ee),'k.'), hold on
        plot(real(active_nodes),imag(active_nodes),'ro')

        plot(z,'gx'), shg
        xlabel('     ')
        axis([-50,15,-30,30])
        title(['sGMRES quadrature rule (k=' num2str(k) ')'], 'FontWeight','normal')
        legend('eigenvalues of A',['generalized evs (m=' num2str(num_it) ')'],['quad nodes (ell=' num2str(ell) ')'])

        % save plot (Figure 5.2 - bottom right)
        mypdf(['fig/network1a_trunc' num2str(k)],.66,1.5)
    end

    %% evaluate quadrature rule for sFOM and plot errors
    %  (note that we are using the GMRES approximate generalized evs for the contour)
    err_sfom_quad = sfom_exp_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_expm,num_it,c,z);
    figure(k-1)
    hold on
    hh1 = line_fewer_markers(1:num_it,err_sfom_quad,10,'s--','Color',[0.8500 0.3250 0.0980]);

    %% sFOM (closed formula)
    err_sfom_closed = sfom_closed_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_expm,@(X) expm(X), num_it);
    figure(k-1)
    h2 = line_fewer_markers(1:num_it,err_sfom_closed,10,'o:','Color',[0.9290 0.6940 0.1250]);

    %% evaluate quadrature rule for sGMRES and plot errors
    err_sgmres = sgmres_exp_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_expm,num_it,c,z);
    figure(k-1)
    h3 = line_fewer_markers(1:num_it,err_sgmres,10,'x-.','Color',[0.4940 0.1840 0.5560]);

    %% funm_quad

    % parameter setup
    param.function = 'exp'; param.restart_length = k;
    param.max_restarts = ceil(num_it/param.restart_length);
    param.tol = 1e-12; param.hermitian = 0; param.V_full = 0;
    param.H_full = 0; param.exact = ex_expm; param.stopping_accuracy = 1e-14;
    param.inner_product = @(a,b) b'*a; param.thick = [];
    param.min_decay = inf; param.waitbar = 1; param.reorth_number = 0;
    param.truncation_length = inf; param.verbose = 1;

    % run restarted Arnoldi
    [f,out] = funm_quad(A,v,param); % number of mat-vec: restart_length*max_restarts

    h4 = line_fewer_markers(param.restart_length*(1:length(out.err)),out.err,10,'--+','Color',[0.4660 0.6740 0.1880]);

    % plot title, axis labels etc.
    title(['convergence (truncation/restart length = ' num2str(k) ')'],'FontWeight','normal')
    xlabel('matrix-vector products m')
    ylabel('2-norm error')
    xlim([0,num_it+1]); ylim([1e-11,1e6]), shg
    legend('best approx','sFOM (quad)','sFOM (closed)','sGMRES','funm\_quad','Location','southwest')
    legend([h0,hh1,h2,h3,h4],'best approx','sFOM (quad)','sFOM (closed)','sGMRES','funm\_quad','Location','SouthWest')

    % save plot (Figure 5.2)
    % k = 2: top left
    % k = 3: top right
    % k = 4: bottom left
    mypdf(['fig/network1b_trunc' num2str(k)],.66,1.5)
    hold off
end