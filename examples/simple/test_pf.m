%% test_pf.m
%  Simple test script for polynomial-filter subroutine.
%  Compute all positive eigenvalues of a random matrix.
%
%  Usage:
%    >> test_pf
%

clear;
N = 2000;
p = 50;

% generate random rank-p matrix
fprintf('Generating matrix (n = %d, p = %d) ...\n', N, p);
X = randn(N, p);
d = randn(1, p);
np = sum(d > 0);

[X, ~] = qr(X, 0);
A = X * diag(d) * X';
A = (A + A') / 2;

%U = randn(N, p);
t1 = tic;
% exact method
opts.maxit = 1000;
opts.tol = 1e-6;
[U, D] = eigs(A, np, 'LA', opts);
D = diag(D);
t1 = toc(t1);

% polynomial-filtered method
t2 = tic;
% initialize PF workspace
work = struct();
work.lb = min(d);
work.ub = 0.2 * work.lb;
work.degree = 8;
S0 = randn(N, ceil(1.1 * np));
work.S0 = S0;
work.ncv = ceil(1.1 * np);

fprintf('   method |    time       | ||AU - UD||_F\n');
fprintf('-------------------------------------------\n');
fprintf('    exact |  %f sec | %e\n', t1, norm(A*U - U * diag(D), 'fro'));

for k = 1:9
    t = tic;
    [U1, D1, work] = eig_inexact(A, work);
    t = toc(t);
    fprintf('   PF(#%d) |  %f sec | %e\n', k, t, norm(A*U1(:, end-np:end) - U1(:, end-np:end)*diag(D1(end-np:end)), 'fro'));
end



