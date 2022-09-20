function y = mylinsolve(M,v,ldl_flag)
% MYLINSOLVE    Solves linear system of equations.
%
% y = mylinsolve(M,v) solves the system M*v = y. The LU
% factors of the coefficient matrix M are reused when the
% same input matrix is identified by its (1,1) element.
%
% To delete the LU factors call MYLINSOLVE with empty argument
% list.


persistent sol indm

verbose = 1;

if nargin < 2
    sol = {}; indm = [];
    if verbose, disp('mylinsolve reset'); end
    return
end
if nargin < 3
    ldl_flag = 0;
end

indm(1) = NaN; % dummy element
m = M(1,1);
ind = find(indm == m,1,'first');

if isempty(ind)
    ind = length(indm) + 1;
    if verbose, disp(['mylinsolve: factorization nr ' num2str(ind-1)]); end
    if ldl_flag
        sol{ind} = decomposition(M,'ldl','upper','CheckCondition',false);
    else
        sol{ind} = decomposition(M);
    end
    indm(ind) = m;
end
if verbose, fprintf(['-' num2str(ind-1)]); end
y = sol{ind}\v;

end
