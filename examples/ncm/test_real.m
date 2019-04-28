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

%% test_real.m
% This is the example of PFPG on the nearest correlation matrix problem
% (NCE). G is a real invalid correlation matrix. See data/*.mat
%
% File obj_comat_v2.m and obj_comat_neg.m illustrate how to insert the
% polynomial-filtered algorithm into the gradient evaluations. Note: the
% negative part of `cor3120.mat` is low-rank.
%
% Usage:
%   >> test_real
%
% See also:
%   obj_comat_v2.m
%   obj_comat_exact.m
%   obj_comat_neg.m
%   obj_comat_exact_neg.m
%   obj_comat_neg_h.m
%%
clear;
addpath data/

opts2 = struct;
opts2.record = 1;
opts2.gtol = 1e-6;

mat_names = {'cor1399', 'cor3120', 'bccd16_t'};
nsize = numel(mat_names);
results_grad = zeros(nsize, 4);
results_pf = zeros(nsize, 4);
results_newton = zeros(nsize, 4);

do_newton = 1;
if ~exist('CorrelationMatrix.m', 'file')
    warning('CorrelationMatrix.m not found. Skipping Newton solvers.');
    do_newton = 0;
end

for ni = 1:nsize
    s = mat_names{ni};
    load(s);
    G = treshape(x, 1);
    G = G + triu(G, 1)';
    n = size(G, 1);
    % opts.record = 1;
    y0 = ones(n, 1) - diag(G);
    [~, g0] = obj_comat_exact(y0, G);
    nrmG0 = norm(g0);
    
    [opts1, obj_exact_h, obj_h] = choose_parameter_real(s);

    tic;
    [y1, ~, out1] = fminGBB(y0, obj_exact_h, opts1, G);
    t1 = toc;
    [f1, g1] = obj_comat_exact(y1, G);
    results_grad(ni, :) = [out1.itr, f1, norm(g1) / nrmG0, t1];
    
    tic;
    [y2, ~, out2] = fminGBB(y0, obj_h, opts1, G);
    t2 = toc;
    [f2, g2] = obj_comat_exact(y2, G);
    results_pf(ni, :) = [out2.itr, f2, norm(g2) / nrmG0, t2];

    if do_newton
        tic;
        [X, y5, out5] = CorrelationMatrix(G, ones(n, 1), 1e-6, 1e-7);
        t5 = toc;
        [f5, g5] = obj_comat_exact(y5, G);
        results_newton(ni, :) = [out5.k, f5, norm(g5) / nrmG0, t5];
    end

end

%% print the result table
print_header('Gradient')
print_result(mat_names, results_grad);
print_header('PFPG');
print_result(mat_names, results_pf);
print_header('Newton');
print_result(mat_names, results_newton);

%% plot the elapsed time
TTable = [results_grad(:, 4), results_pf(:, 4), results_newton(:, 4)];
bar(TTable);
set(gca, 'XTickLabel', mat_names, 'TickLabelInterpreter', 'none');
xlabel('matrix name');
ylabel('time(sec)');
legend('Grad', 'PFPG', 'Newton', 'Location', 'NorthWest');

%% auxiliary functions
function print_result(ns, result)
    fprintf('  name   |  iter      fval      ||g||      time\n');
    fprintf('-----------------------------------------------\n');
    for i = 1:numel(ns)
        fprintf('%8s | %5d   %8.3e  %8.1e  %6.1f \n', ns{i}, result(i,:));
    end
end

function print_header(header)
    fprintf('\n**********************************\n');
    fprintf('* Results of %s\n', header);
    fprintf('**********************************\n');
end

function [opts1, h_exact, h] = choose_parameter_real(mat)
    opts1 = struct;
    opts1.maxit = 1000;
    opts1.tau = 5e-2;
    opts1.ftol = eps;
    opts1.taumax = 2;
    opts1.record = 1;
    opts1.maxcolU = 10;
    opts1.grow = 'none';
    opts1.extrap = 1;
    opts1.exact_cnt = 3;
    opts1.gtol = 1e-5;
    if strcmp(mat, 'cor1399')
        h_exact = @obj_comat_exact;
        h = @obj_comat_v2;
    elseif strcmp(mat, 'cor3120')
        h_exact = @obj_comat_exact;
        h = @obj_comat_neg;
    elseif strcmp(mat, 'bccd16_t')
        opts1.maxcolU = 5;
        h_exact = @obj_comat_exact_neg;
        h = @obj_comat_neg_h;
    else
        error(['cannot determine the function handle for matrix ', mat]);
    end
end
