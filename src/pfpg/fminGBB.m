%
%    PFOpt
%    Copyright (C) 2019  Haoyang Liu (liuhaoyang@pku.edu.cn)
%                        Zaiwen Wen  (wenzw@pku.edu.cn)
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [x, g, out]= fminGBB(x, fun, opts, varargin)
%-------------------------------------------------------------------------
% Line search algorithm for optimization on manifold:
%
%   min f(X), where X \in R^{n,p}
%
%
% Input:
%           X --- 
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= OptManiMulitBallGBB(X0, @fun, opts, data1, data2);
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 mxitr       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
%   
% Output:
%           x --- solution
%           g --- gradient of x
%         Out --- output information
%
%
% Author: Zaiwen Wen
%   Version 1.0 .... 2010/10
%-------------------------------------------------------------------------

% termination rule
if ~isfield(opts, 'xtol');      opts.xtol = 1e-12; end
if ~isfield(opts, 'gtol');      opts.gtol = 1e-5; end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-12; end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'tau');       opts.tau  = 1e-3; end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.1; end
if ~isfield(opts, 'gamma');     opts.gamma  = 0.85; end
if ~isfield(opts, 'STPEPS');    opts.STPEPS  = 1e-10; end
if ~isfield(opts, 'nt');        opts.nt  = 5; end
if ~isfield(opts, 'maxit');     opts.maxit  = 1000; end
if ~isfield(opts, 'eps');       opts.eps = 1e-6; end
if ~isfield(opts, 'record');    opts.record = 0; end
% extrapolations
if ~isfield(opts, 'extrap');    opts.extrap = 1; end
if ~isfield(opts, 'lambda');    opts.lambda = 1e-14; end
if ~isfield(opts, 'maxcolU');   opts.maxcolU = 10; end
% for debugging purpose
if ~isfield(opts, 'debug');     opts.debug = 0; end
% for eig_inexact
if ~isfield(opts, 'k');         opts.k = -1; end
if ~isfield(opts, 'exact_cnt'); opts.exact_cnt = 5; end

%-------------------------------------------------------------------------------
% copy parameters
xtol = opts.xtol;
ftol = opts.ftol;
gtol = opts.gtol;
maxit = opts.maxit;
rhols = opts.rhols;
eta   = opts.eta;
gamma = opts.gamma;
record = opts.record;
nt = opts.nt;
crit = ones(nt, 3);
taumax = opts.taumax;
extrap = opts.extrap;
lambda = opts.lambda;
maxcolU = opts.maxcolU;

[n,p] = size(x);

%% Initial function value and gradient
% prepare for iterations
workproj = struct();
workproj.S0 = [];
workproj.itr = 1;
workproj.force_exact = opts.exact_cnt;
workproj.k = opts.k;
workproj.grow = opts.grow;
workproj.delta = varargin{end};
[f,  g, workproj] = feval(fun, x , varargin{:}, workproj);
out.nfe = 1;  out.fval0 = f;
out.tt = tic;
nrmG  = norm(g, 'fro');
nrmG0 = nrmG;
U = [];
X = [];
if opts.debug; out.err_proj = []; end

Q = 1; Cval = f; tau = opts.tau;
%% Print iteration header if debug == 1
if (record >= 1)
    fprintf('----------- fminBB ----------- \n');
    fprintf('%4s %8s %12s %8s %8s %8s %2s %s\n', ...
        'Iter', 'tau', 'f(X)', 'nrmG', 'xDiff', 'fDiff', 'ls', 'rank');
end

if record == 10; out.fvec = f; end
out.msg = 'exceed max iteration';

%% main iteration
for itr = 1 : maxit
    xp = x;     fp = f;     gp = g;   
    nls = 1; deriv = rhols*nrmG^2;
    
    while 1
        % calculate g, f,
        x = xp - tau*gp;
        %f = feval(fun, x, varargin{:});   out.nfe = out.nfe + 1;
        [f,g,workproj] = feval(fun, x, varargin{:},workproj);   out.nfe = out.nfe + 1;
        if opts.debug
            [~, ~, ww] = obj_comat_exact(x, varargin{:});
            err = norm(ww.B - workproj.B, 2) / norm(ww.B, 2);
            out.err_proj = [out.err_proj, err];
        end
        
        if f <= Cval - tau*deriv || nls >= 5
            %if nls >= 5; workproj.force_exact = 2; end
            break
        end
        tau = eta*tau;
        nls = nls+1;
    end  
    
    % perform extrapolation
    if extrap && itr > 10
        if size(U, 2) < maxcolU
            U = [U, x - xp];
            X = [X, xp];
        else
            c = extrapolation(U, lambda);
            x = X * c;
            [f, g, workproj] = feval(fun, x, varargin{:}, workproj);
            U = [];
            X = [];
        end
    end
    
    if record == 10; out.fvec = [out.fvec; f]; end
    
    nrmG  = norm(g, 'fro');
    s = x - xp; XDiff = norm(s,'fro')/sqrt(n);
    FDiff = abs(fp-f)/(abs(fp)+1);


    if (record >= 1)
        fprintf('%4d %8.2e %12.6e %8.2e %8.2e %8.2e %2d %d/%d\n', ...
            itr, tau, f, nrmG / nrmG0, XDiff, FDiff, nls, workproj.nev, workproj.ncv);
    end
    
    crit(itr,:) = [nrmG, XDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    if nrmG / nrmG0 < 1e-3
        workproj.degree = 8;
    else
        workproj.degree = 5;
    end
    
    %if nrmG < gtol
    if ( FDiff < ftol ) || nrmG / nrmG0 < gtol ...
           % || all(mcrit(2:3) < 10*[xtol, ftol])  
        out.msg = 'converge';
        break;
    end
    
    y = g - gp;      
    sy = abs(iprod(x,y));    tau = opts.tau;
    if sy > 0
        ss = iprod(s, s);
        yy = iprod(y, y);
        
        if mod(itr,2) == 0
            tau = iprod(s, s) / sy;
        else
            tau = sy / iprod(y, y);
        end
        tau = max(ss / sy, sy / yy);
        %tau = iprod(s, s) / sy;
        % safeguarding on tau
        tau = max(min(tau, taumax), 1e-10);
    end
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + f)/Q;
end

out.tt = toc(out.tt);
out.nrmG = nrmG;
out.fval = f;
out.itr = itr;
out.iter = itr;
out.nev = workproj.nev;
end

function a = iprod(x,y)
    a = real(sum(sum(x.*y)));
end

function c = extrapolation(U, lambda)
    k = size(U, 2);
    M = U' * U + lambda * eye(k);

    % solve linear system
    c = M \ ones(k, 1);
    c = c / sum(c);

end
