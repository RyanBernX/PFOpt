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

function [F, g, work_out] = obj_comat_neg_h(y, G, work_in)
    n = numel(y);
    % B = G + Ay
    B = G + spdiags(y, 0, n, n);
    % project B onto S_n^+
    % note: compute the other end
    %work_in.degree = 8;
    work_in.expand_factor = 1.1;
    work_in.shift_factor = 0.01;
    work_in.k = -1;
    work_in.itr = 1;
    work_in.trim = 1;
    work_in.ntrim = 4;
    work_in.h_exact = @eigs_psd_wrapper;
    work_in.h_exact_args = {10, 'LA', 'SubspaceDimension', 40};
    [V, d, work_out] = eig_inexact(-B, work_in);
    id = d > 0;
    B = B + V(:, id) * diag(d(id)) * V(:, id)';
    
    work_out.B = B;
    F = sum(sum(B .* B)) / 2 - sum(y);
    g = diag(B) - 1;
end
