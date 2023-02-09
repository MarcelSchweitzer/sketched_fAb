%% Experiment from "Section 5.2 - Network example"
%
% This script is not optimized for best performance, as its primary
% purpose is generating plots that illustrate the convergence behavior
% of methods.

% Add paths etc.
addpath(genpath('algorithms'),'auxiliary','funm_quad');
mydefaults

% load matrix and solution
load('data/gnutella_comp.mat'); 

% initializations
A = L;
N = size(A,1); 


% Number of Arnoldi iterations to perform
num_it = 150;
[~,~,~,Vfull] = bta(A,speye(N),b,inf(1,num_it),num_it,@(X) X, true, 1);
err_best = best_approx_eval_error(Vfull,ex_gnutella,num_it);

%% truncated Arnoldi

%% truncated Arnoldi
for k = 2:4 % how many blocks to orthogonalise against
    figure(k-1)
    h0 = semilogy(err_best,'k');
    hold on
%     num_it = 120;
    hS = setup_sketching_handle(N,2*num_it); % s = 2*m_max
    [SV,SAV,~,Vtrunc] = bta(A,speye(N),b,inf(1,num_it),k,hS); % number of mat-vec products = num_it
    
    % whitening the basis
    [SV, SAV, Rw] = whiten_basis(SV, SAV);
        
    %% sFOM (closed formula)
    err_sfom_closed = sfom_closed_eval_error(Vtrunc,SV,SAV,hS(b),Rw,ex_gnutella,@(X) sqrtm(X), num_it);
    figure(k-1)
    h2 = line_fewer_markers(1:num_it,err_sfom_closed,10,'o:','Color',[0.9290 0.6940 0.1250]);

    %% adaptive quadrature rule
    [c, z, ell] = invsqrt_adaptive_quad(SV(:,1:num_it),SAV(:,1:num_it),hS(b),1e-10);
    err_sgmres = sgmres_sqrt_quad_eval_error(Vtrunc,SV,SAV,hS(b),Rw,ex_gnutella,num_it,c,z);
    h3 = line_fewer_markers(1:num_it,err_sgmres,23,'x-.','Color',[0.4940 0.1840 0.5560]);

      %% funm_quad

%     % parameter setup
%     param.function = 'exp'; param.restart_length = k;
%     param.max_restarts = ceil(num_it/param.restart_length);
%     param.tol = 1e-12; param.hermitian = 0; param.V_full = 0;
%     param.H_full = 0; param.exact = ex_expm; param.stopping_accuracy = 1e-14;
%     param.inner_product = @(a,b) b'*a; param.thick = [];
%     param.min_decay = inf; param.waitbar = 1; param.reorth_number = 0;
%     param.truncation_length = inf; param.verbose = 1;
% 
%     % run restarted Arnoldi
%     [f,out] = funm_quad(A,v,param); % number of mat-vec: restart_length*max_restarts
% 
%     h4 = line_fewer_markers(param.restart_length*(1:length(out.err)),out.err,10,'--+','Color',[0.4660 0.6740 0.1880]);

    cond_Vtrunc = zeros(1,num_it);
    for i = 1:num_it
        cond_Vtrunc(i) = cond(Vtrunc(:,1:i));
    end
    h5 = semilogy(1:num_it,eps*cond_Vtrunc,'.','Color',[0.3010, 0.7450, 0.9330]);

    % plot title, axis labels etc.
    title(['convergence (truncation length = ' num2str(k) ')'],'FontWeight','normal')
    xlabel('matrix-vector products m')
    ylabel('2-norm error')
    ylim([1e-10,1e1]), shg
    legend([h0,h2,h3,h5],'best approx','sFOM (closed)','sGMRES','u\cdot\kappa(V_m)','Location','SouthWest')
    mypdf(['fig/gnutella' num2str(k)],.66,1.5)
    hold off
end