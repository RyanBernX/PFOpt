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

function [U, d, work_out] = eig_inexact(A, work_in, n)
    %% [U, d, work_out] = eig_inexact(A, work_in)
    %  Perform inexact eigenvalue decomposition of matrix A, using the
    %  the polynomial-filtered algorithm (PF). The output
    %  at least contains the (approximate) positive part or k (approximate)
    %  largest eigen-pairs.
    %  Note: this function is suitable for matrices whos positive part
    %        is low-rank or k is small enough.
    %  Input parameters:
    %    A - input matrix or function handle A*X.
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
    %      work_in.trim - (parameter, experimental) 0 or 1, whether to trim
    %      the other end of the spectra. This option is useful when the
    %      width of the unwanted spectra is very large.
    %
    %      work_in.ntrim - (parameter, experimental) when trim is set to
    %      1, ntrim specifies the number of eigenvalues to be trimmed.
    %
    %      work_in.h_exact - (parameter) users are allowed to provide their
    %      routines for exact eigenvalue decompositions. h_exact should
    %      have the form [U, d] = h_exact(A, k, args), where d is a vector.
    %      The parameters should at least contain A and k (number of
    %      eigenvalues to be computed). Other arguments `args` can be
    %      passed with a cell array, see work_in.h_exact_args.
    %      work_in.h_exact_args - (parameter) if work_in.h_exact is
    %      defined, `h_exact_args` are passed as `args` for `h_exact`.
    %    n - when A is a function handle, n specifies the dimension of A.
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
    if isnumeric(A)
        Afun = @(x) A * x;
        n = size(A, 1);
    else
        Afun = A;
    end
    
    init;
    % initialize work_out
    work_out = work_in;
    work_out.cnt = work_out.cnt + 1;
    

    % check whether we need exact eig
    if work_out.force_exact > 0 || isempty(work_out.S0)
        work_out.force_exact = work_out.force_exact - 1;
        if work_in.k > 0 % for k largest eigenpairs
            if isa(work_in.h_exact, 'function_handle')
                [U, d] = work_in.h_exact(A, work_in.k, work_in.h_exact_args{:});
            else
                opts_eigs = struct();
                opts_eigs.p = min(max(20, 4 * work_out.k), n);
                opts_eigs.tol = 1e-8;
                opts_eigs.issym = 1;
                opts_eigs.isreal = 1;
                %opts_eigs.fail = 'keep';
                [U, D] = eigs(Afun, n, work_in.k, 'LA', opts_eigs);
                d = diag(D); d = d(end:-1:1); U = U(:, end:-1:1);
            end
            % export workspace data
            work_out.ub = d(1);
            work_out.nev = work_in.k;
            work_out.ncv = work_in.k;
            work_out.S0 = U;
        else % for psd and [-1, 1]
            if isa(work_in.h_exact, 'function_handle')
                [U, d] = work_in.h_exact(A, work_in.h_exact_args{:});
            elseif isnumeric(A) && ~issparse(A)
                [U, d] = eig(A, 'vector');
            else
                full_A = Afun(eye(n));
                [U, d] = eig(full_A, 'vector');
            end
            % export workspace data
            if work_in.k == -1
                ln = min(d);
                id = d > 0;
                work_out.lb = 1.1 * ln;
                work_out.ub = work_in.shift_factor * ln;
            elseif work_in.k == -2
                work_out.lb = -0.8;
                work_out.ub = 0.8;
                id = abs(d) > 1;
            end
            work_out.nev = sum(id);
            work_out.ncv = min(size(U, 2), 5 + ceil(work_in.expand_factor * work_out.nev));
            work_out.S0 = U(:, (end - work_out.ncv + 1):end);
        end
        return;
    end

    % check whether to update lb
    if ~work_out.trim && (work_out.cnt >= work_out.nlb || isnan(work_out.lb))
        work_out.cnt = 0;
        % call eigs
        opts_eigs = struct;
        opts_eigs.p = 8;
        opts_eigs.tol = 1e-4;
        opts_eigs.issym = 1;
        opts_eigs.isreal = 0;
        opts_eigs.fail = 'keep';
        if ~isempty(work_out.v0); opts_eigs.v0 = work_out.v0; end
        [work_out.v0, ln] = eigs(Afun, n, 1, 'SR', opts_eigs);
        work_out.lb = 1.1 * ln;
        work_out.shift = work_in.shift_factor * ln;
    end
    
    % check whether we need trimming
    if work_out.trim
        opts_eigs = struct();
        opts_eigs.p = work_in.ntrim * 2;
        opts_eigs.fail = 'keep';
        opts_eigs.tol = 1e-6;
        [Ut, dt] = eigs(Afun, n, work_in.ntrim, 'SR', opts_eigs);
        dt = diag(dt) + 2;
        Afun = @(X) Afun(X) - Ut * (dt .* (Ut' * X));
        work_out.lb = dt(end);
        work_out.ub = work_out.lb * work_in.shift_factor;
    end

    % perform chebyshev iteration
    U = work_out.S0;
    for i=1:work_out.itr
        if isempty(work_out.pcoeff)
            % use cheby
            U = MCheby(Afun, U, work_out.lb, work_out.ub, work_out.degree);
        else
            % use pre-defined poly
            U = polyvalAX(work_out.pcoeff, Afun, U, work_out.lb, work_out.ub);
        end
        % normalize
        colnorm = sqrt(sum(U .^ 2, 1));
        U = U ./ colnorm;
        WW = U' * U;
        %fprintf('%e\n', rcond(WW));
        if rcond(WW) < 1e-6; break; end
    end

    % QR
    if work_out.chol
        % chol-qr
        U = chol_qr(U);
    else
        [U, ~] = qr(U, 0);
    end

    % RR
    AU = Afun(U);
    W = U' * (AU);
    W = (W + W') / 2;

    [V, d] = eig(W, 'vector');
    U = U * V;

    % update work_out
    if work_in.k > 0 % for k largest
        work_out.ub = 0.9 * real(d(1));
        work_out.nev = work_in.k;
        if strcmp(work_in.grow, 'maxeig')
            % multiplicy of lambda_1?
            l1 = max(real(d));
            r1 = sum(d > l1 - 1e-4);
            work_out.k = 10;
        end
    elseif work_in.k == -1 % for psd
        work_out.nev = sum(d > 0);
    elseif work_in.k == -2 % for [-1, 1]
        work_out.nev = sum(abs(d) > 1);
    end
    work_out.ncv = min(n, 5 + ceil(work_in.expand_factor * work_out.nev));

    % check whether need padding
    n_padding = work_in.ncv - work_out.ncv;
    if n_padding < 0
        work_out.S0 = [randn(n, -n_padding), U];
    else
        if work_in.k == -2
            [~, p] = sort(abs(d));
        else
            p = 1:work_in.ncv;
        end
        p = p((1 + n_padding):work_in.ncv);    
        work_out.S0 = U(:, p);
    end

    function init
        if ~isfield(work_in, 'cnt'); work_in.cnt = 0; end
        if ~isfield(work_in, 'force_exact'); work_in.force_exact = 0; end
        if ~isfield(work_in, 'S0'); work_in.S0 = []; end
        if ~isfield(work_in, 'v0'); work_in.v0 = []; end
        if ~isfield(work_in, 'itr'); work_in.itr = 1; end
        if ~isfield(work_in, 'chol'); work_in.chol = 0; end
        if ~isfield(work_in, 'shift'); work_in.shift = 0; end
        if ~isfield(work_in, 'lb'); work_in.lb = nan; end
        if ~isfield(work_in, 'ub'); work_in.ub = nan; end
        if ~isfield(work_in, 'nlb'); work_in.nlb = 8; end
        if ~isfield(work_in, 'degree'); work_in.degree = 8; end
        if ~isfield(work_in, 'pcoeff'); work_in.pcoeff = []; end
        if ~isfield(work_in, 'expand_factor'); work_in.expand_factor = 1.1; end
        if ~isfield(work_in, 'shift_factor'); work_in.shift_factor = 0.2; end
        if ~isfield(work_in, 'k'); work_in.k = -1; end
        if ~isfield(work_in, 'grow'); work_in.grow = 'none'; end
        if ~isfield(work_in, 'nev_max'); work_in.nev_max = floor(0.075 * n); end
        if ~isfield(work_in, 'delta'); work_in.delta = 1e-2; end
        if ~isfield(work_in, 'trim'); work_in.trim = 0; end
        if ~isfield(work_in, 'ntrim'); work_in.ntrim = 5; end
        if ~isfield(work_in, 'h_exact'); work_in.h_exact = 0; end
        if ~isfield(work_in, 'h_exact_args'); work_in.h_exact_args = {}; end
    end
end

function [Q, R] = chol_qr(X)
    H = X' * X;
    R = chol(H);
    % R'R = X'X

    Q = (R' \ X')';
end
