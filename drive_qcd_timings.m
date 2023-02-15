%% Experiment from "Section 5.3 - Lattice QCD example"
%
% This script measures the run time of the different algorithms
% and creates Table 5.1 from the manuscript

% Add paths etc.
addpath(genpath('algorithms'),'auxiliary','funm_quad');
rng(42)

% load matrix and solution
load('data/qcd_matrix_nonhermitian.mat'); load('data/qcd_nonhermitian_exact.mat');

% initializations
A = @(v) Q*(Q*v);
N = size(Q,1); 
b = zeros(N,1); b(1) = 1;
v = Q*b; normv = norm(v); v = v/norm(v);
ex_qcd = 1/normv*qcd_nonhermitian_exact;

num_it = 300;
k = 2;

tol_stop = 1e-5;
tol_quad = tol_stop;

fprintf('Running standard FOM...')
tic
[x_arnoldi, rel_err_arnoldi, ~, iter_arnoldi] = arnoldi_invsqrtm(A, v, num_it, tol_stop, ex_qcd);
time_arnoldi = toc;
fprintf(' done.\n')

fprintf('Running sketched FOM (closed) ...')
tic
[x_fom_closed, rel_err_fom_closed, ~, ~, iter_fom_closed] = sketched_fom_invsqrtm(A,v,num_it,2*num_it+1,k,tol_stop,tol_quad,false,ex_qcd);
time_fom_closed = toc;
fprintf(' done.\n')

fprintf('Running sketched FOM (quad) ...')
tic
[x_fom_quad, rel_err_fom_quad, ~, ~, iter_fom_quad] = sketched_fom_invsqrtm(A,v,num_it,2*num_it+1,k,tol_stop,tol_quad,true,ex_qcd);
time_fom_quad = toc;
fprintf(' done.\n')

fprintf('Running sketched GMRES (quad) ...')
tic
[x_fom_gmres, rel_err_gmres, ~, ~, iter_gmres] = sketched_gmres_invsqrtm(A,v,num_it,2*num_it+1,k,tol_stop,tol_quad,ex_qcd);
time_gmres = toc;
fprintf(' done.\n')

param.function = 'invSqrt';
param.function = 'invSqrt'; param.restart_length = 20;
param.max_restarts = 20;
param.tol = 1e-1*tol_quad; param.hermitian = 0; param.V_full = 0;
param.H_full = 0; param.exact = ex_qcd; param.stopping_accuracy = tol_stop;
param.inner_product = @(a,b) b'*a; param.thick = [];
param.min_decay = inf; param.waitbar = 1; param.reorth_number = 0;
param.truncation_length = inf; param.verbose = 1;
param.transformation_parameter = 1;

fprintf('Running restarted Arnoldi ...')
tic
[x_funm_quad,out] = funm_quad(A,v,param);
time_funmquad = toc;
rel_err_funmquad = norm(x_funm_quad-ex_qcd)/norm(ex_qcd);
iter_funmquad = param.restart_length*length(out.update);
fprintf(' done.\n')

%% print table with results
fprintf('\n\n');
fprintf('           method           | k resp. r | Krylov dim. |  time | rel. err \n')
fprintf('------------------------------------------------------------------------\n')
fprintf('sketched FOM (closed form)  |     %1d     |     %3d     | %1.2fs | %.2e \n', k, iter_fom_closed, time_fom_closed, rel_err_fom_closed);
fprintf('sketched FOM (quadrature)   |     %1d     |     %3d     | %1.2fs | %.2e \n', k, iter_fom_quad, time_fom_quad, rel_err_fom_quad);
fprintf('sketched GMRES (quadrature) |     %1d     |     %3d     | %1.2fs | %.2e \n', k, iter_gmres, time_gmres, rel_err_gmres);
fprintf('funm_quad                   |    %2d     |     %3d     | %1.2fs | %.2e \n', param.restart_length, iter_funmquad, time_funmquad, rel_err_funmquad);
fprintf('standard FOM                |    --     |     %3d     | %1.2fs | %.2e \n\n', iter_arnoldi, time_arnoldi, rel_err_arnoldi);

%% Create LaTeX code for table (Table 5.1)
fID = fopen('fig/tab_qcd.tex','w');
fprintf(fID,'\\begin{table}\n');
fprintf(fID,'\\centering\n');
fprintf(fID,'\\caption{Lattice QCD example. Truncation (resp.\\ restart) length, required number of Krylov iterations, wall-clock time and relative error norm at the final iteration for the different discussed algorithms when invoked with a target accuracy of $10^{-5}$. Details on the experimental setup are given in the final paragraphs of section~\\ref{subsec:qcd}.}\n');
fprintf(fID,'\\label{tab:qcd}\n');
fprintf(fID,'\\begin{tabular}{l|cccc}\n');
fprintf(fID,'method & $k$ resp.\\ $r$ & Krylov dim.\\ & time & rel.\\ error \\\\[1mm]\n');
fprintf(fID,'\\hline\\hline\\\\[-2mm]\n');
fprintf(fID,'sketched FOM (closed form) & %1d & %3d & %1.2fs & $%.2f \\cdot 10^{%1d}$ \\\\\n', k, iter_fom_closed, time_fom_closed, rel_err_fom_closed*10^(-floor(log10(rel_err_fom_closed))),floor(log10(rel_err_fom_closed)));
fprintf(fID,'sketched FOM (quadrature) & %1d & %3d & %1.2fs & $%.2f \\cdot 10^{%1d}$ \\\\\n', k, iter_fom_quad, time_fom_quad, rel_err_fom_quad*10^(-floor(log10(rel_err_fom_quad))),floor(log10(rel_err_fom_quad)));
fprintf(fID,'sketched GMRES (quadrature) & %1d & %3d & %1.2fs & $%.2f \\cdot 10^{%1d}$ \\\\\n', k, iter_gmres, time_gmres, rel_err_gmres*10^(-floor(log10(rel_err_gmres))),floor(log10(rel_err_gmres)));
fprintf(fID,'\\texttt{funm\\_quad} & %2d & %3d & %1.2fs & $%.2f \\cdot 10^{%1d}$ \\\\\n', k, iter_funmquad, time_funmquad, rel_err_funmquad*10^(-floor(log10(rel_err_funmquad))),floor(log10(rel_err_funmquad)));
fprintf(fID,'Standard FOM & -- & %3d & %1.2fs & $%.2f \\cdot 10^{%1d}$ \\\\\n', iter_arnoldi, time_arnoldi, rel_err_arnoldi*10^(-floor(log10(rel_err_arnoldi))),floor(log10(rel_err_arnoldi)));
fprintf(fID,'\\end{tabular}\n');
fprintf(fID,'\\end{table}');
fclose(fID);