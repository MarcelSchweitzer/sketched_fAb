%% Experiment from "Section 5.1 - Convection-diffusion example"
%
% This script is not optimized for best performance, as its primary
% purpose is generating plots that illustrate the convergence behavior
% of methods. 


% Add paths etc.
addpath(genpath('algorithms'),'auxiliary');
mydefaults

% load matrix and solution
load('data/convdiff_matrix.mat'); load('data/convdiff_sol.mat');

% initializations
N = size(A,1); 
v = ones(N,1); v = v/norm(v);

% Number of Arnoldi iterations to perform
num_it = 200; 

%% full Arnoldi and best approx
[~,~,~,Vfull] = bta(A,speye(N),v,inf(1,num_it),num_it,@(X) X);
err_best = best_approx_eval_error(Vfull,ex_convdiff,num_it);
figure(1)
h0 = semilogy(err_best,'k');
hold on


%% truncated Arnoldi
k = 4; % how many blocks to orthogonalise against
hS = setup_sketching_handle(N,2*num_it); % s = 2*m_max
[SV,SAV,~,Vtrunc] = bta(A,speye(N),v,inf(1,num_it),k,hS); % number of mat-vec products = num_it

% whitening the basis
[SV, SAV, Rw] = whiten_basis(SV, SAV);

%% adaptive quadrature rule
%  one to rule them all (determined for m_max)
[c, z, ell] = invsqrt_adaptive_quad(SV(:,1:num_it),SAV(:,1:num_it),hS(v),1e-7);

%% evaluate quadrature rule for sFOM and plot errors
%  (note that we are using the GMRES approximate generalized evs for the contour)
err_sfom_quad = sfom_invsqrt_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_convdiff,num_it,c,z);
figure(1)
hh1 = line_fewer_markers(1:num_it,err_sfom_quad,23,'s--','Color',[0.8500 0.3250 0.0980]);


%% sFOM (closed formula)
[err_sfom_closed, diff_f] = sfom_closed_eval_error_and_diff(Vtrunc,SV,SAV,A,Vfull,hS(v),Rw,ex_convdiff,@(X) inv(sqrtm(full(X))), num_it);
figure(1)
h2 = line_fewer_markers(1:num_it,err_sfom_closed,23,'o:','Color',[0.9290 0.6940 0.1250]);

% for plotting Ritz values in order to better understand irregular
% convergence behavior
ritz_fom84 = eig((SV(:,1:84)'*SAV(:,1:84))/(SV(:,1:84)'*SV(:,1:84)));
ritz_fom85 = eig((SV(:,1:85)'*SAV(:,1:85))/(SV(:,1:85)'*SV(:,1:85)));

%% evaluate quadrature rule for sGMRES and plot errors
err_sgmres = sgmres_invsqrt_quad_eval_error(Vtrunc,SV,SAV,hS(v),Rw,ex_convdiff,num_it,c,z);
figure(1)
h3 = line_fewer_markers(1:num_it,err_sgmres,23,'x-.','Color',[0.4940 0.1840 0.5560]);

% plot vertical line between 84th and 85th iteration
h4 = semilogy([84.5,84.5],[1e-16,1e16],'-','Color',[.5,.5,.5]);
uistack(h4,'bottom')

% plot conditioning of truncated basis
cond_Vtrunc = zeros(1,num_it);
for i = 1:num_it
    cond_Vtrunc(i) = cond(Vtrunc(:,1:i));
end
h5 = semilogy(1:num_it,eps*cond_Vtrunc,'.','Color',[0.4660, 0.6740, 0.1880]);

h6 = semilogy(1:num_it,diff_f,'-.','Color',[0.3010, 0.7450, 0.9330]);

% plot title, axis labels etc.
title(['convergence (truncation length k=' num2str(k) ')'],'FontWeight','normal')
xlabel('matrix-vector products m')
ylabel('2-norm error')
ylim([1e-16,1e3]), shg
xlim([0,num_it+1])
legend([h0,hh1,h2,h3,h5],'best approx','sFOM (quad)','sFOM (closed)','sGMRES','u\cdot\kappa(V_m)','Location','SouthWest')

% save plot (Figure 5.1 - left)
mypdf(['fig/convdiff_trunc' num2str(k)],.66,1.5)

%% plot Ritz values, eigenvalues and quadrature nodes
figure
eA = eigs(A,20,'smallestabs');
hhhh1 = plot(eA,'k.');
hold on
hhhh2 = plot(ritz_fom84,'b+');
hhhh3 = plot(ritz_fom85,'ro');
nodes = -(1-z)./(1+z);
hhhh4 = plot(real(nodes),imag(nodes),'gx');
axis([-10,50,-120,120]), shg
% replot the red dots on top to make fully visible
plot(ritz_fom85,'ro');

title('sketched FOM Ritz values (k=4)','FontWeight','normal')
xlabel('   ')
legend([hhhh1,hhhh2,hhhh3,hhhh4],'eigenvalues of A','Ritz values (m=84)','Ritz values (m=85)','quadrature nodes','Location','SouthWest')

% save plot (Figure 5.1 - right)
mypdf(['fig/convdiffritz_trunc' num2str(k)],.66,1.5)