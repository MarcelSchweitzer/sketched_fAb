function [SV,SAV,SBV,Vfull] = bta(A,B,v,xi,t,hS,mgs,reo)
% block truncated rational Arnoldi with randomized reduction
% (A,B) - matrix pair
% v     - block vector
% xi    - shift parameters (can be infinite)
% t     - truncation parameter (t=2 will give Lanczos for Hermitian AB)
% hS    - handle to left-side basis reduction
% 
% SV    - left-reduced Krylov basis matrix
% SAV   - left-reduced A times basis matrix
% SBV   - left-reduced B times basis matrix
% Vfull - also return full unreduced Krylov basis (optional)

[N,b] = size(v);

if nargin < 7
    mgs = false;
end
if nargin < 8
    reo = 0;
end
if ishermitian(A) && ishermitian(B) && all(isreal(xi))
    disp('Hermitian problem with real shifts. (LDL flag = 1, only for mylinsolve.m)')
    ldl_flag = 1;
else
    ldl_flag = 0;
end

if norm(B - speye(N),'fro') < eps
    disp('B is identity matrix.')
    identityB = 1;
else
    identityB = 0;
end


V = zeros(N, t*b); % truncated orthonormal basis (last t blocks)

mylinsolve(); % using decomposition, slower but more accurate
%mylinsolve_umfpack(); % fast, less accuracte, used by eigs?
%util_linsolve(); % LU, faster but less accurate
[v,~] = qr(v,0);
V = [ V(:,b+1:end) , v ];
SV = hS(v); s = size(SV,1);
SV = [ SV, zeros(s,length(xi)*b) ];
SAV = zeros(s,(length(xi)+1)*b); 
SBV = zeros(s,(length(xi)+1)*b);
Vfull = [];
if nargout>=4
    Vfull = zeros(N,(length(xi)+1)*b);
    Vfull(:,1:b) = v;
end
wb = waitbar(0,'bta running');
for j = 1:length(xi)
    w = V(:,end-b+1:end);

    Aw = A*w; 

    % compute these retrospectively
    SAV(:,1+(j-1)*b:(j)*b) = hS(Aw); 
    if identityB
        Bw = w;
        SBV(:,1+(j-1)*b:(j)*b) = SV(:,1+(j-1)*b:(j)*b); 
    else
        Bw = B*w;
        SBV(:,1+(j-1)*b:(j)*b) = hS(Bw); 
    end

    if isfinite(xi(j))
        w = mylinsolve(A - xi(j)*B, Bw, ldl_flag); % slower
        %w = mylinsolve_umfpack(A - xi(j)*B, Bw); % used by eigs?
        %w = util_linsolve(A - xi(j)*B, Bw); 
    else
        if identityB
            w = Aw;
        else
            w = mylinsolve(B, Aw, ldl_flag);
            %w = mylinsolve_umfpack(B, Aw); % used by eigs?
            %w = util_linsolve(B, Aw);
        end
    end
    if mgs % modified gram schmidt
        for reo = 0:1
            for i = 1:t
                w = w - V(:,1+(i-1)*b:i*b)*(V(:,1+(i-1)*b:i*b)'*w);
            end
        end
    else % classical gram schmidt
        w = w - V*(V'*w);
    end
    [v,~] = qr(w,0);
    V = [ V(:,b+1:end) , v ];
    
    SV(:,1+j*b:(j+1)*b) = hS(v); 
    
    if nargout>=4
        Vfull(:,1+j*b:(j+1)*b) = v;
    end
    waitbar(j/length(xi),wb)
end

% final block column of SAV and SBV
SAV(:,1+j*b:(j+1)*b) = hS(A*v); 
if identityB
    SBV(:,1+j*b:(j+1)*b) = SV(:,1+j*b:(j+1)*b);
else
    SBV(:,1+j*b:(j+1)*b) = hS(B*v);
end

close(wb)
end

