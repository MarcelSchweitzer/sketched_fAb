%% Experiment from "Section 5.3 - Lattice QCD example"
%
% This script is not optimized for best performance, as its primary
% purpose is generating plots that illustrate the convergence behavior
% of methods.

% Add paths etc.
addpath(genpath('algorithms'),'auxiliary','funm_quad');
mydefaults

% load matrix and solution
load('data/qcd_matrix_nonhermitian.mat'); load('data/qcd_nonhermitian_exact.mat');

% initializations
A = Q*Q; % explicitly forming Q^2 is normally not advisable of course!
N = size(A,1); 
b = zeros(N,1); b(1) = 1;
v = Q*b; normv = norm(v); v = v/norm(v);
ex_qcd = 1/normv*qcd_nonhermitian_exact;

% Number of Arnoldi iterations to perform
num_it = 300;

%% full Arnoldi and best approx
[~,~,~,Vfull] = bta(A,speye(N),v,inf(1,num_it),num_it,@(X) X);
err_best = best_approx_eval_error(Vfull,ex_qcd,num_it);

%% truncated Arnoldi
for k = 2:4 % how many blocks to orthogonalise against
    figure(k-1)
    h0 = semilogy(err_best,'k');
    hold on

    hS = setup_sketching_handle(N,2*num_it); % s = 2*m_max
    [SV,SAV,~,Vtrunc] = bta(A,speye(N),v,inf(1,num_it),k,hS); % number of mat-vec products = num_it

    % whitening the basis
    [SV, SAV, Rw] = whiten_basis(SV, SAV);

    %% adaptive quadrature rule
    %  one to rule them all (determined for m_max)
    [c, z, ell] = invsqrt_adaptive_quad(SV(:,1:num_it),SAV(:,1:num_it),hS(v),1e-10);

    %% evaluate quadrature rule for sFOM and plot errors
    %  (note that we are using the GMRES approximate generalized evs for the contour)
    err_sfom_quad = sfom_invsqrt_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_qcd,num_it,c,z);
    figure(k-1)
    hold on
    hh1 = line_fewer_markers(1:num_it,err_sfom_quad,23,'s--','Color',[0.8500 0.3250 0.0980]);

    %% sFOM (closed formula)
    err_sfom_closed = sfom_closed_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_qcd,@(X) inv(sqrtm(full(X))), num_it);
    figure(k-1)
    h2 = line_fewer_markers(1:num_it,err_sfom_closed,23,'o:','Color',[0.9290 0.6940 0.1250]);

    %% evaluate quadrature rule for sGMRES and plot errors
    err_sgmres = sgmres_invsqrt_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_qcd,num_it,c,z);
    figure(k-1)
    h3 = line_fewer_markers(1:num_it,err_sgmres,23,'x-.','Color',[0.4940 0.1840 0.5560]);

    %% funm_quad

    % parameter setup
    param.function = 'invSqrt'; param.restart_length = k;
    param.max_restarts = ceil(num_it/param.restart_length);
    param.tol = 1e-13; param.hermitian = 0; param.V_full = 0;
    param.H_full = 0; param.exact = ex_qcd; param.stopping_accuracy = 1e-14;
    param.inner_product = @(a,b) b'*a; param.thick = [];
    param.min_decay = inf; param.waitbar = 1; param.reorth_number = 0;
    param.truncation_length = inf; param.verbose = 1;
    param.transformation_parameter = 1;

    % run restarted Arnoldi
    [f,out] = funm_quad(A,v,param); % number of mat-vec: restart_length*max_restarts

    h4 = line_fewer_markers(param.restart_length*(1:length(out.err)),out.err,10,'--+','Color',[0.4660 0.6740 0.1880]);

    % plot title, axis labels etc.
    title(['convergence (truncation/restart length = ' num2str(k) ')'],'FontWeight','normal')
    xlabel('matrix-vector products m')
    ylabel('2-norm error')
    ylim([1e-7,2e0]), shg
    xlim([0,num_it+1])
    legend([h0,hh1,h2,h3,h4],'best approx','sFOM (quad)','sFOM (closed)','sGMRES','funm\_quad','Location','SouthWest')

    % save plot (Figure 5.3)
    % k = 2: top left
    % k = 3: top right
    % k = 4: bottom left
    mypdf(['fig/qcd_trunc' num2str(k)],.66,1.5)
    hold off

    % also generate a plot with restart length 20
    if k == 2
        param.restart_length = 20;
        param.max_restarts = ceil(num_it/param.restart_length);

        % run restarted Arnoldi
        [f,out] = funm_quad(A,v,param); % number of mat-vec: restart_length*max_restarts

        figure(4)
        h0 = semilogy(err_best,'k');
        hold on
        hh1 = line_fewer_markers(1:num_it,err_sfom_quad,23,'s--','Color',[0.8500 0.3250 0.0980]);
        h2 = line_fewer_markers(1:num_it,err_sfom_closed,23,'o:','Color',[0.9290 0.6940 0.1250]);
        h3 = line_fewer_markers(1:num_it,err_sgmres,23,'x-.','Color',[0.4940 0.1840 0.5560]);
        h4 = line_fewer_markers(param.restart_length*(1:length(out.err)),out.err,10,'--+','Color',[0.4660 0.6740 0.1880]);

        title('convergence (trunc. len = 2 / rest. len. = 20)','FontWeight','normal')
        xlabel('matrix-vector products m')
        ylabel('2-norm error')
        ylim([1e-7,2e0]), shg
        legend([h0,hh1,h2,h3,h4],'best approx','sFOM (quad)','sFOM (closed)','sGMRES','funm\_quad','Location','SouthWest')
        
        % save plot (Figure 5.3 - bottom right)
        mypdf('fig/qcd_trunc2_restart20',.66,1.5)
        hold off
    end
end