function [x,flag,relres,iter,resvec,resveccg] = minres(A,b,tol,maxit,M1,M2,x0,mask,varargin)
%MINRES   Minimum Residual Method.
%   X = MINRES(A,B) attempts to find a minimum norm residual solution X to
%   the system of linear equations A*X=B. The N-by-N coefficient matrix A
%   must be symmetric but need not be positive definite. The right hand
%   side column vector B must have length N.
%
%   X = MINRES(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = MINRES(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then MINRES uses the default, 1e-6.
%
%   X = MINRES(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then MINRES uses the default, min(N,20).
%
%   X = MINRES(A,B,TOL,MAXIT,M) and X = MINRES(A,B,TOL,MAXIT,M1,M2) use
%   symmetric positive definite preconditioner M or M=M1*M2 and effectively
%   solve the system inv(sqrt(M))*A*inv(sqrt(M))*Y = inv(sqrt(M))*B for Y
%   and then return X = inv(sqrt(M))*Y. If M is [] then a preconditioner is
%   not applied.  M may be a function handle returning M\X.
%
%   X = MINRES(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0
%   is [] then MINRES uses the default, an all zero vector.
%
%   [X,FLAG] = MINRES(A,B,...) also returns a convergence FLAG:
%    0 MINRES converged to the desired tolerance TOL within MAXIT
%    iterations.
%    1 MINRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 MINRES stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during MINRES became
%      too small or too large to continue computing.
%    5 preconditioner M was not symmetric positive definite.
%
%   [X,FLAG,RELRES] = MINRES(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = MINRES(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = MINRES(A,B,...) also returns a vector of
%   estimates of the MINRES residual norms at each iteration, including
%   NORM(B-A*X0).
%
%   [X,FLAG,RELRES,ITER,RESVEC,RESVECCG] = MINRES(A,B,...) also returns a
%   a vector of estimates of the Conjugate Gradients residual norms at each
%   iteration.
%
%   Example:
%      n = 100; on = ones(n,1); A = spdiags([-2*on 4*on -2*on],-1:1,n,n);
%      b = sum(A,2); tol = 1e-10; maxit = 50; M = spdiags(4*on,0,n,n);
%      x = minres(A,b,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-------------------------------%
%      function y = afun(x,n)
%      y = 4 * x;
%      y(2:n) = y(2:n) - 2 * x(1:n-1);
%      y(1:n-1) = y(1:n-1) - 2 * x(2:n);
%      %-------------------------------%
%   as input to MINRES:
%      x1 = minres(@(x)afun(x,n),b,tol,maxit,M);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, PCG, QMR, SYMMLQ,
%   TFQMR, ICHOL, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2012 The MathWorks, Inc.
%   $Revision: 1.6.4.13 $ $Date: 2012/02/09 20:58:50 $

if (nargin < 2)
    error(message('MATLAB:minres:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [ma,n] = size(A);
    if (ma ~= n)
        error(message('MATLAB:minres:NonSquareMatrix'));
    end
    if ~isequal(size(b),[ma,1])
        error(message('MATLAB:minres:RSHsizeMatchCoeffMatrix', ma));
    end
else
    ma = size(b,1);
    n = ma;
    if ~iscolumn(b)
        error(message('MATLAB:minres:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol <= eps
    warning(message('MATLAB:minres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:minres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if nargout >= 6, resveccg = 0; end;	% resveccg(1) = norm(b-A*xcg) = norm(0)
    if (nargout < 2)
        itermsg('minres',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[ma,ma])
            error(message('MATLAB:minres:WrongPrecondSize', ma));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[ma,ma])
            error(message('MATLAB:minres:WrongPrecondSize', ma));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

existM = existM1 | existM2;

if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:minres:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 8) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:minres:TooManyInputs'));
end

% Set up for the method
flag = 1;
iter = 0;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,mask,varargin{:});
normr = norm(r);                   % Norm of residual

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    resvec = normr;
    if nargout >= 6, resveccg = normr;  end;
    if (nargout < 2)
        itermsg('minres',tol,maxit,0,flag,iter,relres);
    end
    return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for MINRES residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
if nargout >= 6
    resveccg = zeros(maxit+2,1);   % Preallocate vector for CG residuals
    resveccg(1) = normr;           % resveccg(1) = norm(b-A*x0)
end
normrmin = normr;                  % Norm of minimum residual

vold = r;
if existM1
    u = iterapp('mldivide',m1fun,m1type,m1fcnstr,vold,mask,varargin{:});
    if ~all(isfinite(u))
        flag = 2;
        relres = normr / n2b;
        resvec = resvec(1);
        if nargout >= 6, resveccg = resveccg(1);    end
        if nargout < 2
            itermsg('minres',tol,maxit,0,flag,iter,relres);
        end
        return
    end
else % no preconditioner
    u = vold;
end
if existM2
    v = iterapp('mldivide',m2fun,m2type,m2fcnstr,u,mask,varargin{:});
    if ~all(isfinite(v))
        flag = 2;
        relres = normr / n2b;
        resvec = resvec(1);
        if nargout >= 6, resveccg = resveccg(1);    end
        if nargout < 2
            itermsg('minres',tol,maxit,0,flag,iter,relres);
        end
        return
    end
else % no preconditioner
    v = u;
end
beta1 = vold' * v;
if (beta1 <= 0)
    flag = 5;
    relres = normr / n2b;
    resvec = resvec(1);
    if nargout >= 6, resveccg = resveccg(1); end
    if nargout < 2
        itermsg('minres',tol,maxit,0,flag,iter,relres);
    end
    return
end
beta1 = sqrt(beta1);
snprod = beta1;
vv = v / beta1;
v = iterapp('mtimes',afun,atype,afcnstr,vv,mask,varargin{:});
Amvv = v;
alpha = vv' * v;
v = v - (alpha/beta1) * vold;

% Local reorthogonalization
numer = vv' * v;
denom = vv' * vv;
v = v - (numer/denom) * vv;
volder = vold;
vold = v;

if existM1
    u = iterapp('mldivide',m1fun,m1type,m1fcnstr,vold,mask,varargin{:});
    if ~all(isfinite(u))
        flag = 2;
        relres = normr / n2b;
        resvec = resvec(1);
        if nargout >= 6, resveccg = resveccg(1);    end
        if nargout < 2
            itermsg('minres',tol,maxit,0,flag,iter,relres);
        end
        return
    end
else % no preconditioner
    u = vold;
end
if existM2
    v = iterapp('mldivide',m2fun,m2type,m2fcnstr,u,mask,varargin{:});
    if ~all(isfinite(v))
        flag = 2;
        relres = normr / n2b;
        resvec = resvec(1);
        if nargout >= 6, resveccg = resveccg(1);    end
        if nargout < 2
            itermsg('minres',tol,maxit,0,flag,iter,relres);
        end
        return
    end
else % no preconditioner
    v = u;
end
betaold = beta1;
beta = vold' * v;
if (beta < 0)
    flag = 5;
    relres = normr / n2b;
    resvec = resvec(1);
    if nargout >= 6, resveccg = resveccg(1);    end
    if nargout < 2
        itermsg('minres',tol,maxit,0,flag,iter,relres);
    end
    return
end
iter = 1;
beta = sqrt(beta);
gammabar = alpha;
epsilon = 0;
deltabar = beta;
gamma = sqrt(gammabar^2 + beta^2);
mold = zeros(n,1);
Amold = mold;
m = vv / gamma;
Am = Amvv / gamma;
cs = gammabar / gamma;
sn = beta / gamma;
x = x + snprod * cs * m;
snprodold = snprod;
snprod = snprod * sn;

% This recurrence produces CG iterates.
% Enable the following statement to see xcg.
%xcg = x + snprod * (sn/cs) * m;

if existM
    r = r - snprodold * cs * Am;
    normr = norm(r);
else
    normr = abs(snprod);
end
resvec(2,1) = normr;
if nargout >= 6
    if (cs == 0)
        % It's possible that this cs value is zero (CG iterate does not exist)
        normrcg = Inf;
    else
        rcg = r - snprod*(sn/cs)*Am;
        normrcg = norm(rcg);
    end
    resveccg(2,1) = normrcg;
end

% Check for convergence after first step.
if normr <= tolb
    flag = 0;
    relres = normr / n2b;
    resvec = resvec(1:2);
    if nargout >= 6,    resveccg = resveccg(1:2);    end
    if (nargout < 2)
        itermsg('minres',tol,maxit,1,flag,iter,relres);
    end
    return
end

stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 2 : maxit
    
    vv = v * (1/beta);
    v = iterapp('mtimes',afun,atype,afcnstr,vv,mask,varargin{:});
    Amolder = Amold;
    Amold = Am;
    Am = v;
    v = v - (beta / betaold) * volder;
    alpha = vv' * v;
    v = v - (alpha / beta) * vold;
    volder = vold;
    vold = v;
    if existM1
        u = iterapp('mldivide',m1fun,m1type,m1fcnstr,vold,mask,varargin{:});
        if ~all(isfinite(u))
            flag = 2;
            break
        end
    else % no preconditioner
        u = vold;
    end
    if existM2
        v = iterapp('mldivide',m2fun,m2type,m2fcnstr,u,mask,varargin{:});
        if ~all(isfinite(v))
            flag = 2;
            break
        end
    else % no preconditioner
        v = u;
    end
    betaold = beta;
    beta = vold' * v;
    if (beta < 0)
        flag = 5;
        break
    end
    beta = sqrt(beta);
    delta = cs * deltabar + sn * alpha;
    molder = mold;
    mold = m;
    m = vv - delta * mold - epsilon * molder;
    Am = Am - delta * Amold - epsilon * Amolder;
    gammabar = sn * deltabar - cs * alpha;
    epsilon = sn * beta;
    deltabar = - cs * beta;
    gamma = sqrt(gammabar^2 + beta^2);
    m = m / gamma;
    Am = Am / gamma;
    cs = gammabar / gamma;
    sn = beta / gamma;
    % Check for stagnation of the method
    if (snprod*cs == 0) || (abs(snprod*cs)*norm(m) < eps*norm(x))
        % increment the number of consecutive iterates which are the same
        stag = stag + 1;
    else
        stag = 0;
    end
    x = x + (snprod * cs) * m;
    snprodold = snprod;
    snprod = snprod * sn;
    % This recurrence produces CG iterates.
    % Enable the following statement to see xcg.
    %xcg = x + snprod * (sn/cs) * m;
    
    if existM
        r = r - snprodold*cs*Am;
        normr = norm(r);
    else
        normr = abs(snprod);
    end
    resvec(ii+1,1) = normr;
    if nargout >= 6
        % It's possible that this cs value is zero (CG iterate does not exist).
        if (cs == 0)
            normrcg = Inf;
        else
            rcg = r - snprod*(sn/cs)*Am;
            normrcg = norm(rcg);
        end
        resveccg(ii+2,1) = normrcg;
    end
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        % double check residual norm is less than tolerance
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,mask,varargin{:});
        normr = norm(r);
        resvec(ii+1,1) = normr;
        if (normr <= tolb)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:minres:tooSmallTolerance'));
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    
    if (normr < normrmin)      % update minimal norm quantities
        normrmin = normr;
        xmin = x;
        imin = ii;
    end
    
    if (stag >= maxstagsteps)      % 3 iterates are the same
        flag = 3;
        break
    end
end                                % for ii = 1 : maxit
if isempty(ii)
  ii = 1;
end

% returned solution is first with minimal residual
if (flag == 0)
    relres = normr / n2b;
else
    r_comp = b - iterapp('mtimes',afun,atype,afcnstr,xmin,mask,varargin{:});
    if norm(r_comp) <= normr
        x = xmin;
        iter = imin;
        relres = norm(r_comp) / n2b;
    else
        iter = ii;
        relres = normr / n2b;
    end
end

% truncate the zeros from resvec
if ((flag <= 1) || (flag == 3))
    resvec = resvec(1:ii+1);
    if nargout >= 6,    resveccg = resveccg(1:ii+2);    end
else
    resvec = resvec(1:ii);
    if nargout >= 6,    resveccg = resveccg(1:ii+1);      end
end

% only display a message if the output flag is not used
if (nargout < 2)
    itermsg('minres',tol,maxit,ii,flag,iter,relres);
end
