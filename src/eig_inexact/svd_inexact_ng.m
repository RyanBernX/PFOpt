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

function [U, d, V, work_out] = svd_inexact_ng(Afun, work_in, m, n)
    %% [U, d, V, work_out] = svd_inexact(Afun, work_in, m, n)
    %  Perform inexact eigenvalue decomposition of matrix A, the output
    %  at least contains the (approximate) positive part or k (approximate)
    %  largest eigen-pairs.
    %  Note: this function is suitable for matrices whos positive part
    %        is low-rank or k is small enough.
    %  Input parameters:
    %    A - input matrix.
    %    work_in - workspace of this function. In iterative methods, the
    %      information of previous update can serve as a good warm start of the
    %      current update. Therefore the workspace can be passed along these
    %      updates to obtain better performances. `work_in` is a struct that
    %      contains parameters and data.
    %
    %      work_in.k - (parameter) number of desired eigenvalues to be
    %      computed. Usually k largest. If k == -1 then eig_inexact
    %      computes all positive part.
    %
    %      work_in.cnt - (data) number of current calls. Used for determining
    %      whether to update lb (see work_in.nlb below). Users are not
    %      suggested to modify this field.
    %
    %      work_in.nlb - (parameter) number of calls before update lb using
    %      the smallest eigenvalue of A. The accuracy of lb has great impact
    %      on the effectiveness of polynomial filtering. Default is 8.
    %
    %      work_in.force_exact - (data) force this subroutine to use exact
    %      eigenvalue decomposition [U D] = eig(A). After `eig` is called,
    %      lb, S0, nev, ncv will be set accordingly. Users can set this field
    %      to any positive integer. This field will decrease by 1 on exit.
    %      For example, when work_in.force_exact = 2, this subroutine will use
    %      exact eigenvalue decomposition in the following 2 consecutive calls.
    %      When set to inf, eig_inexact will behave just like eig.
    %
    %      work_in.S0 - (data) initial subspace. It can either be extracted
    %      from exact eigenvalue decomposition or come from the previous
    %      update. Users are not suggested to modify this field.
    %
    %      work_in.nev - (data) the number of positive eigenvalues in S0.
    %      work_in.ncv - (data) the dimension of S0.
    %
    %      work_in.v0 - (data) initial vector for calling `eigs`.
    %
    %      work_in.itr - (parameter) the number of Chebyshev iterations when
    %      performing subspace update. Default is 1. Note: the subspace update
    %      can terminate earlier than normal when then condition number of U'U
    %      is too large.
    %
    %      work_in.chol - (parameter) whether to perform Cholesky-QR when
    %      extracting orthogonal basis from U. Default is 0. Note: Cholesky-QR
    %      is usually faster than QR but less accuate.
    %
    %      work_in.shift - (data) shift relative to zero. Must be non-positive.
    %      work_in.lb - (data) the lower bound of all eigenvalues of A. Must
    %      be negative. If work_in.lb is larger than min(eig(A)), then this
    %      subroutine will fail. If work_in.lb is MUCH SMALLER than min(eig(A))
    %      then the accuracy will be poor.
    %      work_in.ub - (data) the upper bound of all eigenvalues of A to
    %      be suppressed. Automatically set in all scenarios.
    %
    %      work_in.degree - (parameter) the degree of the Chebyshev polynomial.
    %      Ignored when work_in.pcoeff is not empty. Default is 8.
    %
    %      work_in.pcoeff - (parameter) the coefficients of the generic
    %      polynomial. When given, this subroutine will use this polynomial
    %      instead of Chebyshev polynomials.
    %
    %  Output parameters:
    %    U - orthogonal matrix n by work_out.ncv, approximate eigenvectors.
    %    D - vector of size work_out.ncv, approximate eigenvalues.
    %        Note: usually D contains all positive eigenvalues of A.
    %    work_out - updated workspace. The parameters are unchanged. Possible
    %      changed fields are:
    %      work_out.S0 - part of U or augmented from U.
    %      work_out.nev - number of positive eigenvalues of S0/U.
    %      work_out.ncv - dimension of work_out.S0.
    %      work_out.cnt - usually increased by 1 or zeroed out (eigs called).
    %      work_out.v0 - updated when eigs called.
    %      work_out.force_exact - usually decreased by 1.
    %%
    % A can be an operator
    if isnumeric(Afun)
        A = @(x) Afun * x;
        At = @(x) Afun' * x;
        Afun = @(x, s) Afun_proto(x, s, A, At);
    else
        A = @(x) Afun(x, 'notransp');
        At = @(x) Afun(x, 'transp');
    end
    
    init;
    % initialize work_out
    work_out = work_in;
    work_out.cnt = work_out.cnt + 1;
    

    % check whether we need exact eig
    if work_out.force_exact > 0 || isempty(work_out.S0)
        work_out.force_exact = work_out.force_exact - 1;
        if work_in.k > 0 % for k largest eigenpairs
            opts_svds = struct();
            opts_svds.p = min(max(20, 4 * work_out.k), min(m, n));
            opts_svds.tol = 1e-8;
            opts_svds.fail = 'keep';
            [U, D, V] = svds(Afun, [m n], work_in.k, 'largest', opts_svds);
            d = diag(D);
            % export workspace data
            work_out.ub = d(work_in.k);
            work_out.nev = work_in.k;
            work_out.ncv = work_in.k;
            work_out.S0 = U;
            return;
        end
    end
    % check whether we need to expand S0
    n_padding = work_in.k - size(work_in.S0, 2);
    if n_padding > 0
        work_out.S0 = [randn(m, n_padding), work_out.S0];
    end

    % perform chebyshev iteration
    U = work_out.S0;
    AAt = @(x) A(At(x));
    for i=1:work_out.itr
        if isempty(work_out.pcoeff)
            % use cheby
            U = MCheby(AAt, U, work_out.lb, work_out.ub, work_out.degree);
        else
            % use pre-defined poly
            U = polyvalAX(work_out.pcoeff, AAt, U, work_out.lb, work_out.ub);
        end
        % normalize
        colnorm = sqrt(sum(U .^ 2, 1));
        U = U ./ colnorm;
        WW = U' * U;
        %fprintf('%e\n', rcond(WW));
        if rcond(WW) < 1e-6; break; end
    end
    % retrieve V
    % Note: At * p(AAt)*U0 = p(AtA)*(AtU0)
    V = At(U);

    % QR
    if work_out.chol
        % chol-qr
        U = chol_qr(U);
        V = chol_qr(V);
    else
        [U, ~] = qr(U, 0);
        [V, ~] = qr(V, 0);
    end

    % RR
    AV = A(V);
    W = U' * (AV);
    [Y, S, Z] = svd(W, 'econ');
    
    % assemble Ritz pairs
    U = U * Y;
    V = V * Z;
    d = diag(S);

    % update work_out
    if work_in.k > 0 % for k largest
        work_out.ub = d(work_in.k);
        work_out.nev = work_in.k;
    end
    work_out.ncv = min([m, n, 5 + ceil(work_in.expand_factor * work_out.nev)]);

    % check whether need padding
    n_padding = work_in.ncv - work_out.ncv;
    if n_padding < 0
        work_out.S0 = [randn(m, -n_padding), U];
    else
        work_out.S0 = U(:, 1:work_out.ncv);
    end

    function init
        if ~isfield(work_in, 'cnt'); work_in.cnt = 0; end
        if ~isfield(work_in, 'force_exact'); work_in.force_exact = 0; end
        if ~isfield(work_in, 'S0'); work_in.S0 = []; end
        if ~isfield(work_in, 'v0'); work_in.v0 = []; end
        if ~isfield(work_in, 'itr'); work_in.itr = 1; end
        if ~isfield(work_in, 'chol'); work_in.chol = 0; end
        if ~isfield(work_in, 'shift'); work_in.shift = 0; end
        if ~isfield(work_in, 'lb'); work_in.lb = 0; end
        if ~isfield(work_in, 'ub'); work_in.ub = nan; end
        if ~isfield(work_in, 'nlb'); work_in.nlb = 8; end
        if ~isfield(work_in, 'degree'); work_in.degree = 8; end
        if ~isfield(work_in, 'pcoeff'); work_in.pcoeff = []; end
        if ~isfield(work_in, 'expand_factor'); work_in.expand_factor = 1.1; end
        if ~isfield(work_in, 'shift_factor'); work_in.shift_factor = 0.2; end
        if ~isfield(work_in, 'k'); work_in.k = 6; end
    end
end

function [Q, R] = chol_qr(X)
    H = X' * X;
    R = chol(H);
    % R'R = X'X

    Q = (R' \ X')';
end

function Y = Afun_proto(X, trans, A, At)
    if strcmp(trans, 'transp')
        Y = At(X);
    else
        Y = A(X);
    end
end
