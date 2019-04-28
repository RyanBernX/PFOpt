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

%% test_phaselift.m
% Example codes for the phase recovery problem, using the phase lift
% formulation and the GAUGE solver. See https://epubs.siam.org/doi/10.1137/15M1034283
%
% The demo codes compare the performance of GAUGE (by Michael P. Friedlander)
% and PFGAUGE (polynomial-filtered GAUGE). The PFGAUGE implementation is
% based on GAUGE. See http://www.cs.ubc.ca/ ∼ mpf/low-rank-opt
%
% Usage:
%   >> test_phaselift.m
%
% See also:
%   low-rank-opt/saga_sd.m
%   low-rank-opt/+hop/pl.m
%
%%

low_rank_opt_dir = './low-rank-opt';
addpath(low_rank_opt_dir);
run([low_rank_opt_dir '/setpath.m']);

clear;

solverOpts.verbosity  = 0;
solverOpts.iterations = 200;
solverOpts.eigTolMin = 1e-6;
solverOpts.eigInexactDegree = 5;

n = 2048;
iseed = 20;

xerrs = zeros(1, 2);
nffts = zeros(1, 2);
iters = zeros(1, 2);
times = zeros(1, 2);

fprintf('------------------------------------------------------------------------------\n');
fprintf('    |            GAUGE                   |              PFGAUGE               \n');
fprintf('    --------------------------------------------------------------------------\n');
fprintf('  L |  time  | itr |   nFFT   |   xerr   |  time  | itr |   nFFT   |   xerr   \n');
fprintf('------------------------------------------------------------------------------\n');

for L = 12:-1:7
    % GAUGE
    [A, b, x0, info] = ...
        experiments.gendatapl('n',n,'m',1,'L',L,'seed',iseed);

    solverOpts.eigInexact = 0;
    [x, r, info_sol] = saga_sd(A, b, solverOpts);

    if info_sol.stat == 1
        xError = util.hermitianerror(x(:),x0(:),'fro');
        xErrorRel = xError / norm(x0)^2;

        xerrs(1) = xErrorRel;
        nffts(1) = info_sol.nfft;
        iters(1) = info_sol.nit;
        times(1) = info_sol.time;

    else
        xerrs(1) = inf;
        nffts(1) = inf;
        iters(1) = inf;
        times(1) = inf;
    end

    % PFGAUGE
    [A, b, x0, info] = ...
       experiments.gendatapl('n',n,'m',1,'L',L,'seed',iseed);

    solverOpts.eigInexact = 1;
    [x, r, info_sol] = saga_sd(A, b, solverOpts);

    xError = util.hermitianerror(x(:),x0(:),'fro');
    xErrorRel = xError / norm(x0)^2;

    if info_sol.stat == 1
        xError = util.hermitianerror(x(:),x0(:),'fro');
        xErrorRel = xError / norm(x0)^2;

        xerrs(2) = xErrorRel;
        nffts(2) = info_sol.nfft;
        iters(2) = info_sol.nit;
        times(2) = info_sol.time;
    else
        xerrs(2) = inf;
        nffts(2) = inf;
        iters(2) = inf;
        times(2) = inf;
    end
    fprintf('%3d | %6.2f | %3d | %8.1e | %8.1e | %6.2f | %3d | %8.1e | %8.1e \n', ...
        L, times(1), iters(1), nffts(1), xerrs(1), ...
           times(2), iters(2), nffts(2), xerrs(2));
end

