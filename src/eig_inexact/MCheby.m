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

function X = MCheby(Afun,X,a,b,degree)
% X = MCheby(Afun,X,a,b,degree)
%
if degree == 0; return; end
if degree == 1; X = Afun(X); return; end
%
c0 = (a + b)/(a - b);
c1 = 2 / (b - a);
B = @(X) c0*X + c1*Afun(X);
%
T1 = X; T2 = B(X);
for j = 2:degree
    X = 2*B(T2) - T1;
    if j < degree
        T1 = T2; T2 = X;
    end
end
end
