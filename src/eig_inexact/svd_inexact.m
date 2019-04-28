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

function [U, d, V, work_out] = svd_inexact(Afun, work_in, m, n)
    % A can be numeric or handle
    if isnumeric(Afun)
        A = @(x) Afun * x;
        At = @(x) Afun' * x;
    else
        A = @(x) Afun(x, 'notransp');
        At = @(x) Afun(x, 'transp');
    end
    AA = @(x) applyH(A, At, x, m, n);
    [UU, d, work_out] = eig_inexact(AA, work_in, m + n);
    
    % only keep the positive half
    [d, ipiv] = sort(d, 'descend');
    UU = UU(:, ipiv);
    id = d > 0;
    d = d(id);
    U = sqrt(2) * UU(1:m, id);
    V = sqrt(2) * UU(m+1:m+n, id);
end

function Y = applyH(A, At, X, m, n)
    Y = zeros(size(X));
    Y(1:m, :) = A(X(m+1:m+n, :));
    Y(m+1:m+n, :) = At(X(1:m, :));
end
