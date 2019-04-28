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

function [F, g, holder] = obj_comat_exact(y, G, varargin)
    n = numel(y);
    % B = G + Ay
    B = G + spdiags(y, 0, n, n);
    % project B onto S_n^+
    [V, d] = eig(B, 'vector');
    id = d > 0;
    B = V(:, id) * diag(d(id)) * V(:, id)';
    
    F = sum(sum(B .* B)) / 2 - sum(y);
    g = diag(B) - 1;
    holder = struct();
    holder.nev = sum(id);
    holder.ncv = sum(id);
    holder.B = B;
end
