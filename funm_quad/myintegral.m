function [Q,errbnd,npts] = myintegral(fun,a,b,varargin)
%MYINTEGRAL  Numerically evaluate integral.
%
%   NOTE: This is a minor modification of Matlab's INTEGRAL routine, which
%   returns the number of adaptively chosen quadrature nodes.
% 
%   Q = MYINTEGRAL(FUN,A,B) approximates the integral of function FUN from 
%   A to B using global adaptive quadrature and default error tolerances.
%
%   FUN must be a function handle. A and B can be -Inf or Inf. If both are
%   finite, they can be complex. If at least one is complex, MYINTEGRAL
%   approximates the path integral from A to B over a straight line path.
%
%   For scalar-valued problems the function Y = FUN(X) must accept a vector
%   argument X and return a vector result Y, the integrand function
%   evaluated at each element of X. For array-valued problems (see the
%   'ArrayValued' option below) FUN must accept a scalar and return an
%   array of values.
%
%   Q = MYINTEGRAL(FUN,A,B,PARAM1,VAL1,PARAM2,VAL2,...) performs the
%   integration with specified values of optional parameters. The available
%   parameters are
%
%   'AbsTol', absolute error tolerance
%   'RelTol', relative error tolerance
%
%       MYINTEGRAL attempts to satisfy |Q - I| <= max(AbsTol,RelTol*|Q|),
%       where I denotes the exact value of the integral. Usually RelTol
%       determines the accuracy of the integration. However, if |Q| is
%       sufficiently small, AbsTol determines the accuracy of the
%       integration, instead. The default value of AbsTol is 1.e-10, and
%       the default value of RelTol is 1.e-6. Single precision integrations
%       may require larger tolerances.
%
%   'ArrayValued', FUN is an array-valued function when the input is scalar
%
%       When 'ArrayValued' is true, FUN is only called with scalar X, and
%       if FUN returns an array, MYINTEGRAL computes a corresponding array of
%       outputs Q. The default value is false.
%
%   'Waypoints', vector of integration waypoints
%
%       If FUN(X) has discontinuities in the interval of integration, the
%       locations should be supplied as a 'Waypoints' vector. Waypoints
%       should not be used for singularities in FUN(X). Instead, split the
%       interval and add the results from separate integrations with
%       singularities at the endpoints. If A, B, or any entry of the
%       waypoints vector is complex, the integration is performed over a
%       sequence of straight line paths in the complex plane, from A to the
%       first waypoint, from the first waypoint to the second, and so
%       forth, and finally from the last waypoint to B.
%
%   Examples:
%       % Integrate f(x) = exp(-x^2)*log(x)^2 from 0 to infinity:
%       f = @(x) exp(-x.^2).*log(x).^2
%       Q = myintegral(f,0,Inf)
%
%       % To use a parameter in the integrand:
%       f = @(x,c) 1./(x.^3-2*x-c)
%       Q = myintegral(@(x)f(x,5),0,2)
%
%       % Specify tolerances:
%       Q = myintegral(@(x)log(x),0,1,'AbsTol',1e-6,'RelTol',1e-3)
%
%       % Integrate f(z) = 1/(2z-1) in the complex plane over the
%       % triangular path from 0 to 1+1i to 1-1i to 0:
%       Q = myintegral(@(z)1./(2*z-1),0,0,'Waypoints',[1+1i,1-1i])
%
%       % Integrate the vector-valued function sin((1:5)*x) from 0 to 1:
%       Q = myintegral(@(x)sin((1:5)*x),0,1,'ArrayValued',true)
%
%   Class support for inputs A, B, and the output of FUN:
%      float: double, single
%
%   See also INTEGRAL2, INTEGRAL3, FUNCTION_HANDLE

%   Portions based on "quadva" by Lawrence F. Shampine.
%   Ref: L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
%   Journal of Computational and Applied Mathematics 211, 2008, pp.131-140

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/11/13 04:01:00 $
%
%   Modifications to be used in conjunction with the quadrature-based
%   restart code FUNM_QUAD described in 
%
%   A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%   stable Arnoldi restarts for matrix functions based on quadrature,
%   SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
%   Don't use otherwise.



% Validate the first three inputs.
%narginchk(3,inf);
if ~isa(fun,'function_handle')
    error(message('MATLAB:integral:funArgNotHandle'));
end
if ~(isscalar(a) && isfloat(a) && isscalar(b) && isfloat(b))
    error(message('MATLAB:integral:invalidEndpoint'));
end
opstruct = integralParseArgs(varargin{:});
[Q,errbnd,npts] = integralCalc(fun,a,b,opstruct);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = Gauss7Kronrod15
% Gauss-Kronrod (7,15) pair.
nodes = [ ...
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, ...
    -0.7415311855993944, -0.5860872354676911, -0.4058451513773972, ...
    -0.2077849550078985, 0, 0.2077849550078985, ...
    0.4058451513773972, 0.5860872354676911, 0.7415311855993944, ...
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126];
wt15 = [ ...
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502, ...
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854, ...
    0.2044329400752989, 0.2094821410847278, 0.2044329400752989, ...
    0.1903505780647854, 0.1690047266392679, 0.1406532597155259, ...
    0.1047900103222502, 0.06309209262997855, 0.02293532201052922];
wt7 = [ ...
    0, 0.1294849661688697, 0, ...
    0.2797053914892767, 0, 0.3818300505051189, ...
    0, 0.4179591836734694, 0, ...
    0.3818300505051189, 0, 0.2797053914892767, ...
    0, 0.1294849661688697, 0];
r = struct('Nodes',nodes,'HighWeights',wt15,'LowWeights',wt7);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [q,errbnd,npts] = integralCalc(fun,a,b,opstruct)
%INTEGRALCALC  Perform INTEGRAL calculation.
%
%   NOTE: This is a minor modification of Matlab's internal INTEGRALCALC 
%   routine, which returns the number of adaptively chosen quadrature 
%   nodes.
% 

npts = 0;

%   Based on "quadva" by Lawrence F. Shampine.
%   Ref: L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
%   Journal of Computational and Applied Mathematics 211, 2008, pp.131-140
%
%   Strong algebraic singularities are handled with a technique from
%   L.F. Shampine, "Weighted Quadrature by Change of Variable", Neural,
%   Parallel, and Scientific Computations 18, 2010, pp. 195-206.
%
%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/11/13 04:01:02 $

% Variable names in all caps are referenced in nested functions.
% Fixed parameters.
DEFAULT_MAXINTERVALCOUNT = 16384; % 2^14
% Unpackage options.
rule = opstruct.Rule;
NODES = rule.Nodes(:);
NNODES = numel(NODES);
WT = rule.HighWeights;
EWT = WT - rule.LowWeights;
ATOL = opstruct.AbsTol;
RTOL = opstruct.RelTol;
INITIALINTERVALCOUNT = opstruct.InitialIntervalCount;
ARRAYVALUED = opstruct.ArrayValued;
MAXINTERVALCOUNT = ceil(DEFAULT_MAXINTERVALCOUNT*opstruct.Persistence);
waypoints = opstruct.Waypoints(:).';
FUN = fun;
if ~(isreal(a) && isreal(b) && isreal(waypoints))
    % Handle contour integration.
    tinterval = [a,waypoints,b];
    if any(~isfinite(tinterval))
        error(message('MATLAB:integral:nonFiniteContourError'));
    end
    A = a;
    B = b;
    [q,errbnd] = vadapt(@identityTransform,tinterval);
else
    % Define A and B and note the direction of integration on real axis.
    if b < a
        % Integrate left to right and change sign at the end.
        reversedir = true;
        A = b;
        B = a;
    else
        reversedir = false;
        A = a;
        B = b;
    end
    % Trim and sort the waypoints vector.
    waypoints = sort(waypoints(waypoints>A & waypoints<B));
    % Construct interval vector with waypoints.
    interval = [A,waypoints,B];
    % Extract A and B from interval vector to regularize possible mixed
    % single/double inputs.
    A = interval(1);
    B = interval(end);
    ZERO = zeros(class(interval));
    ONE = ones(class(interval));
    % Identify the task and perform the integration.
    if A == B
        % Handles both finite and infinite cases.
        % Return zero or nan of the appropriate class.
        q = midpArea(FUN,A,B);
        errbnd = q;
    elseif isfinite(A) && isfinite(B)
        if numel(interval) > 2
            % Analytical transformation suggested by K.L. Metlov:
            alpha = 2*sin( asin((A + B - 2*interval(2:end-1))/(A - B))/3 );
            interval = [-ONE,alpha,ONE];
        else
            interval = [-ONE,ONE];
        end
        [q,errbnd] = vadapt(@AtoBInvTransform,interval);
    elseif isfinite(A) && isinf(B)
        if numel(interval) > 2
            alpha = sqrt(interval(2:end-1) - A);
            interval = [ZERO,alpha./(ONE+alpha),ONE];
        else
            interval = [ZERO,ONE];
        end
        [q,errbnd] = vadapt(@AToInfInvTransform,interval);
    elseif isinf(A) && isfinite(B)
        if numel(interval) > 2
            alpha = sqrt(B - interval(2:end-1));
            interval = [-ONE,-alpha./(ONE+alpha),ZERO];
        else
            interval = [-ONE,ZERO];
        end
        [q,errbnd] = vadapt(@minusInfToBInvTransform,interval);
    elseif isinf(A) && isinf(B)
        if numel(interval) > 2
            % Analytical transformation suggested by K.L. Metlov:
            alpha = tanh( asinh(2*interval(2:end-1))/2 );
            interval = [-ONE,alpha,ONE];
        else
            interval = [-ONE,ONE];
        end
        if isa(A,'single') || isa(B,'single')
            interval = single(interval);
        end
        [q,errbnd] = vadapt(@minusInfToInfInvTransform,interval);
    else % i.e., if isnan(a) || isnan(b)
        q = midpArea(FUN,A,B);
        errbnd = q;
    end
    % Account for integration direction.
    if reversedir
        q = -q;
    end
end

%--------------------------------------------------------------------------

    function [q,errbnd] = vadapt(u,tinterval)
        % Iterative routine to perform the integration.
        % Compute the path length and split tinterval if needed.
        [tinterval,pathlen] = split(tinterval,INITIALINTERVALCOUNT);
        if pathlen == 0
            % Test case: integral(@(x)x,1+1i,1+1i);
            q = midpArea(FUN,A,B);
            errbnd = q;
            return
        end
        % Set MAXINTERVALCOUNT high enough that each initial subinterval
        % can be refined at least once.
        MAXINTERVALCOUNT = max(MAXINTERVALCOUNT,2*(numel(tinterval)-1));
        if ARRAYVALUED
            [q,errbnd] = iterateArrayValued(u,tinterval,pathlen);
        else
            [q,errbnd] = iterateScalarValued(u,tinterval,pathlen);
        end
    end % vadapt

%--------------------------------------------------------------------------

    function [q,errbnd] = iterateArrayValued(u,tinterval,pathlen)
        firstFunEval = true;
        % Initialize array of subintervals of [a,b].
        subs = [tinterval(1:end-1);tinterval(2:end)];
        nsubs = size(subs,2);
        % Initialize main loop
        while true
            % SUBS contains subintervals of [a,b] where the integral is not
            % sufficiently accurate. The first row of SUBS holds the left
            % end points and the second row, the corresponding right
            % endpoints.
            midpt = sum(subs)/2;   % midpoints of the subintervals
            halfh = diff(subs)/2;  % half the lengths of the subintervals
            x = bsxfun(@plus,NODES*halfh,midpt); % NNODES x nsubs
            if firstFunEval
                % Unwrap the first function evaluation to get size and
                % class information.
                [t,w] = u(x(:,1));
                fxj = FUN(t(1)).*w(1);
                nfx = numel(fxj);
                outsize = size(fxj);
                outcls = superiorfloat(x,fxj);
                finalInputChecks(x,fxj);
                % Define initializers.
                if issparse(fxj)
                    zeros_nfx_by_nsubs = sparse(nfx,nsubs);
                    zeros_outsize = sparse(outsize(1),outsize(2));
                    zeros_nfx_by_1 = sparse(nfx,1);
                else
                    zeros_nfx_by_nsubs = zeros([nfx,nsubs],outcls);
                    zeros_outsize = zeros(outsize,outcls);
                    zeros_nfx_by_1 = zeros(nfx,1,outcls);
                end
                qsubs = zeros_nfx_by_nsubs;
                errsubs = zeros_nfx_by_nsubs;
                qsubsk = fxj*WT(1);
                errsubsk = fxj*EWT(1);
                % Finish up the first subinterval.
                for j = 2:NNODES
                    fxj = FUN(t(j)).*w(j);
                    qsubsk = qsubsk + fxj*WT(j);
                    errsubsk = errsubsk + fxj*EWT(j);
                end
                qsubsk = qsubsk.*halfh(1);
                errsubsk = errsubsk.*halfh(1);
                q = qsubsk;
                qsubs(:,1) = qsubsk(:);
                errsubs(:,1) = errsubsk(:);
                % Now process the remaining subintervals. We don't check
                % the mesh spacing on the first iteration.
                for k = 2:nsubs
                    qsubsk = zeros_outsize;
                    errsubsk = zeros_outsize;
                    [t,w] = u(x(:,k));
                    for j = 1:NNODES
                        fxj = FUN(t(j)).*w(j);
                        qsubsk = qsubsk + fxj*WT(j);
                        errsubsk = errsubsk + fxj*EWT(j);
                    end
                    qsubsk = qsubsk.*halfh(k);
                    errsubsk = errsubsk.*halfh(k);
                    q = q + qsubsk;
                    qsubs(:,k) = qsubsk(:);
                    errsubs(:,k) = errsubsk(:);
                end
                firstFunEval = false;
                % Initialize partial sums.
                q_ok = zeros_nfx_by_1;
                err_ok = zeros_nfx_by_1;
            else
                qsubs = zeros_nfx_by_nsubs;
                errsubs = zeros_nfx_by_nsubs;
                q = reshape(q_ok,outsize);
                for k = 1:nsubs
                    qsubsk = zeros_outsize;
                    errsubsk = zeros_outsize;
                    [t,w] = u(x(:,k));
                    too_close = checkSpacing(t);
                    if too_close
                        break
                    end
                    for j = 1:NNODES
                        fxj = FUN(t(j)).*w(j);
                        qsubsk = qsubsk + fxj*WT(j);
                        errsubsk = errsubsk + fxj*EWT(j);
                    end
                    qsubsk = qsubsk.*halfh(k);
                    errsubsk = errsubsk.*halfh(k);
                    q = q + qsubsk;
                    qsubs(:,k) = qsubsk(:);
                    errsubs(:,k) = errsubsk(:);
                end
                if too_close
                    break
                end
            end
            % Locate subintervals where the approximate integrals are
            % sufficiently accurate and use them to update the partial
            % error sum.
            nleft = 0;
            tol = RTOL*abs(q(:));
            tolr = 2*tol/pathlen;
            tola = 2*ATOL/pathlen;
            err_not_ok = zeros_nfx_by_1;
            for k = 1:nsubs
                abserrsubsk = abs(errsubs(:,k));
                if all(abserrsubsk <= tolr*halfh(k) | ...
                        abserrsubsk <= tola*halfh(k))
                    q_ok = q_ok + qsubs(:,k);
                    err_ok = err_ok + errsubs(:,k);
                else
                    err_not_ok = err_not_ok + abserrsubsk;
                    nleft = nleft + 1;
                    subs(:,nleft) = subs(:,k);
                    midpt(nleft) = midpt(k);
                end
            end
            
            npts = npts + length(WT)*(nsubs-nleft);
            
            subs = subs(:,1:nleft); % Trim subs array.
            
            midpt = midpt(1:nleft); % Trim midpt array.
            % The approximate error bound is constructed by adding the
            % approximate error bounds for the subintervals with accurate
            % approximations to the 1-norm of the approximate error bounds
            % for the remaining subintervals.  This guards against
            % excessive cancellation of the errors of the remaining
            % subintervals.
            errbnd = reshape(abs(err_ok) + err_not_ok, outsize);
            % Check for nonfinites.
            if ~all(isfinite(q(:)) & isfinite(errbnd(:)))
                warning(message('MATLAB:integral:NonFiniteValue'));
                if opstruct.ThrowOnFail
                    error(message('MATLAB:integral:unsuccessful'));
                end
                break
            end
            % Test for convergence.
            if nleft < 1 || all(errbnd(:) <= tol | errbnd(:) <= ATOL)
                break
            end
            % Split the remaining subintervals in half. Quit if splitting
            % results in too many subintervals.
            nsubs = 2*nleft;
            if nsubs > MAXINTERVALCOUNT
                % newPersistence = ceil(nsubs/DEFAULT_MAXINTERVALCOUNT);
                warning(message('MATLAB:integral:MaxIntervalCountReached', ...
                    sprintf('%9.1e',max(errbnd(:)))));
                if opstruct.ThrowOnFail
                    error(message('MATLAB:integral:unsuccessful'));
                end
                break
            end
            subs = reshape([subs(1,:); midpt; midpt; subs(2,:)],2,[]);
        end
    end

%--------------------------------------------------------------------------

    function [q,errbnd] = iterateScalarValued(u,tinterval,pathlen)
        firstFunEval = true;
        % Initialize array of subintervals of [a,b].
        subs = [tinterval(1:end-1);tinterval(2:end)];
        % Initialize partial sums.
        q_ok = 0;
        err_ok = 0;
        % Initialize main loop
        while true
            % SUBS contains subintervals of [a,b] where the integral is not
            % sufficiently accurate. The first row of SUBS holds the left
            % end points and the second row, the corresponding right
            % endpoints.
            midpt = sum(subs)/2;   % midpoints of the subintervals
            halfh = diff(subs)/2;  % half the lengths of the subintervals
            x = bsxfun(@plus,NODES*halfh,midpt); % NNODES x nsubs
            x = reshape(x,1,[]);
            [t,w] = u(x); % Transform back to the original domain.
            if firstFunEval
                fx = FUN(t);
                finalInputChecks(x,fx);
                fx = fx.*w;
                firstFunEval = false;
            else
                too_close = checkSpacing(t);
                if too_close
                    break
                end
                fx = FUN(t).*w;
            end
            fx = reshape(fx,numel(WT),[]);
            % Quantities for subintervals.
            qsubs = (WT*fx) .* halfh;
            errsubs = (EWT*fx) .* halfh;
            % Calculate current values of q and tol.
            q = sum(qsubs) + q_ok;
            tol = max(ATOL,RTOL*abs(q));
            % Locate subintervals where the approximate integrals are
            % sufficiently accurate and use them to update the partial
            % error sum.
            ndx = find(abs(errsubs) <= (2*tol/pathlen)*halfh);
            npts = npts + length(WT)*length(ndx);
            err_ok = err_ok + sum(errsubs(ndx));
            % Remove errsubs entries for subintervals with accurate
            % approximations.
            errsubs(ndx) = [];
            % The approximate error bound is constructed by adding the
            % approximate error bounds for the subintervals with accurate
            % approximations to the 1-norm of the approximate error bounds
            % for the remaining subintervals.  This guards against
            % excessive cancellation of the errors of the remaining
            % subintervals.
            errbnd = abs(err_ok) + norm(errsubs,1);
            % Check for nonfinites.
            if ~(isfinite(q) && isfinite(errbnd))
                warning(message('MATLAB:integral:NonFiniteValue'));
                if opstruct.ThrowOnFail
                    error(message('MATLAB:integral:unsuccessful'));
                end
                break
            end
            % Test for convergence.
            if errbnd <= tol
                break
            end
            % Remove subintervals with accurate approximations.
            subs(:,ndx) = [];
            % Update the partial sum for the integral.
            q_ok = q_ok + sum(qsubs(ndx));
            midpt(ndx) = []; % Remove unneeded midpoints.
            if isempty(subs)
                break
            end
            % Split the remaining subintervals in half. Quit if splitting
            % results in too many subintervals.
            nsubs = 2*size(subs,2);
            if nsubs > MAXINTERVALCOUNT
                % newPersistence = ceil(nsubs/DEFAULT_MAXINTERVALCOUNT);
                warning(message('MATLAB:integral:MaxIntervalCountReached', ...
                    sprintf('%9.1e',errbnd)));
                if opstruct.ThrowOnFail
                    error(message('MATLAB:integral:unsuccessful'));
                end
                break
            end
            subs = reshape([subs(1,:); midpt; midpt; subs(2,:)],2,[]);
        end
    end

%--------------------------------------------------------------------------

    function q = midpArea(f,a,b)
        % Return q = (b-a)*f((a+b)/2). Although formally correct as a low
        % order quadrature formula, this function is only used to return
        % nan or zero of the appropriate class when a == b, isnan(a), or
        % isnan(b).
        x = (a + b)/2;
        if isfinite(a) && isfinite(b) && ~isfinite(x)
            % Treat overflow, e.g. when finite a and b > realmax/2
            x = a/2 + b/2;
        end
        fx = f(x);
        if ~all(isfinite(fx(:)))
            warning(message('MATLAB:integral:NonFiniteValue'));
            if opstruct.ThrowOnFail
                error(message('MATLAB:integral:unsuccessful'));
            end
       end
        q = (b - a)*fx;
    end % midpArea

%--------------------------------------------------------------------------

    function [x,w] = identityTransform(t)
        x = t;
        w = ones(size(t),class(t));
    end

%--------------------------------------------------------------------------

    function [x,w] = AtoBInvTransform(t)
        % Transform to weaken singularities at both ends: [a,b] -> [-1,1]
        x = 0.25*(B-A)*t.*(3 - t.^2) + 0.5*(B+A);
        w = 0.75*(B-A)*(1 - t.^2);
    end

%--------------------------------------------------------------------------

    function [x,w] = AToInfInvTransform(t)
        % Transform to weaken singularity at left end: [a,Inf) -> [0,Inf).
        % Then transform to finite interval: [0,Inf) -> [0,1].
        tt = t ./ (1 - t);
        x = A + tt.^2;
        w = 2*tt ./ (1 - t).^2;
    end

%--------------------------------------------------------------------------

    function [x,w] = minusInfToBInvTransform(t)
        % Transform to weaken singularity at right end: (-Inf,b] -> (-Inf,b].
        % Then transform to finite interval: (-Inf,b] -> (-1,0].
        tt = t ./ (1 + t);
        x = B - tt.^2;
        w = -2*tt ./ (1 + t).^2;
    end

%--------------------------------------------------------------------------

    function [x,w] = minusInfToInfInvTransform(t)
        % Transform to finite interval: (-Inf,Inf) -> (-1,1).
        x = t ./ (1 - t.^2);
        w = (1 + t.^2) ./ (1 - t.^2).^2;
    end

%--------------------------------------------------------------------------

    function too_close = checkSpacing(x)
        ax = abs(x);
        tcidx = find(abs(diff(x)) <= 100*eps(class(x))* ...
            max(ax(1:end-1),ax(2:end)),1);
        too_close = ~isempty(tcidx);
        if too_close
            warning(message('MATLAB:integral:MinStepSize', ...
                num2str(x(tcidx),6)));
            if opstruct.ThrowOnFail
                error(message('MATLAB:integral:unsuccessful'));
            end
        end
    end % checkSpacing

%--------------------------------------------------------------------------

    function [x,pathlen] = split(x,minsubs)
        % Split subintervals in the interval vector X so that, to working
        % precision, no subinterval is longer than 1/MINSUBS times the
        % total path length. Removes subintervals of zero length, except
        % that the resulting X will always has at least two elements on
        % return, i.e., if the total path length is zero, X will be
        % collapsed into a single interval of zero length.  Also returns
        % the integration path length.
        absdx = abs(diff(x));
        if isreal(x)
            pathlen = x(end) - x(1);
        else
            pathlen = sum(absdx);
        end
        if minsubs > 1
            if pathlen > 0
                udelta = minsubs/pathlen;
                nnew = ceil(absdx*udelta) - 1;
                idxnew = find(nnew > 0);
                nnew = nnew(idxnew);
                for j = numel(idxnew):-1:1
                    k = idxnew(j);
                    nnj = nnew(j);
                    % Calculate new points.
                    newpts = x(k) + (1:nnj)./(nnj+1)*(x(k+1)-x(k));
                    % Insert the new points.
                    x = [x(1:k),newpts,x(k+1:end)];
                end
            end
            % Remove useless subintervals.
            x(abs(diff(x))==0) = [];
            if isscalar(x)
                % Return at least two elements.
                x = [x(1),x(1)];
            end
        end
    end % split

%--------------------------------------------------------------------------

    function finalInputChecks(x,fx)
        % Do final input validation with sample input and outputs to the
        % integrand function.
        % Check classes.
        if ~isfloat(fx)
            error(message('MATLAB:integral:UnsupportedClass',class(fx)));
        end
        % Check sizes.
        if ~ARRAYVALUED && ~isequal(size(x),size(fx))
            error(message('MATLAB:integral:FxNotSameSizeAsX'));
        end
        outcls = superiorfloat(x,fx);
        outdbl = strcmp(outcls,'double');
        % Make sure that RTOL >= 100*eps(outcls) except when using pure
        % absolute error control (ATOL>0 && RTOL==0). Obviously we're
        % assuming that the default tolerances are OK.
        if ~(ATOL > 0 && RTOL == 0) && RTOL < 100*eps(outcls)
            RTOL = 100*eps(outcls);
        end
        if outdbl
            % Single RTOL or ATOL should not force any single precision
            % computations.
            RTOL = double(RTOL);
            ATOL = double(ATOL);
        else
            WT = cast(WT,class(fx));
            EWT = cast(EWT,class(fx));
            NODES = cast(NODES,class(x));
        end
    end % finalInputChecks
end % integralCalc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = integralParseArgs(varargin)
%INTEGRALPARSEARGS  Parse optional arguments to INTEGRAL.

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/09/23 19:06:37 $

p = inputParser;
p.addParamValue('AbsTol',1e-10,@validateAbsTol);
p.addParamValue('RelTol',1e-6,@validateRelTol);
p.addParamValue('Waypoints',[],@validateWaypoints);
p.addParamValue('ArrayValued',false,@validateArrayValued);
p.parse(varargin{:});
s = p.Results;
s.Rule = Gauss7Kronrod15;
s.InitialIntervalCount = 10;
s.Persistence = 1;
s.ThrowOnFail = false;
end

%--------------------------------------------------------------------------

function p = validateAbsTol(x)
if ~(isfloat(x) && isscalar(x) && isreal(x) && x >= 0)
    error(message('MATLAB:integral:invalidAbsTol'));
end
p = true;
end

%--------------------------------------------------------------------------

function p = validateRelTol(x)
if ~(isfloat(x) && isscalar(x) && isreal(x) && x >= 0)
    error(message('MATLAB:integral:invalidRelTol'));
end
p = true;
end

%--------------------------------------------------------------------------

function p = validateWaypoints(x)
if ~(isvector(x) || isequal(x,[])) || any(~isfinite(x))
    error(message('MATLAB:integral:invalidWaypoints'));
end
p = true;
end

%--------------------------------------------------------------------------

function p = validateArrayValued(x)
if ~(isscalar(x) && islogical(x))
    error(message('MATLAB:integral:invalidArrayValued'));
end
p = true;

%--------------------------------------------------------------------------
end

