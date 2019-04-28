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
%% install_newton.m
%  Download the most recent CorrelationMatrix.m from
%  http://www.math.nus.edu.sg/~matsundf/CorrelationMatrix.m
%  NOTE: THE SOURCE CODE ONLINE IS MISSING!

url = 'http://www.math.nus.edu.sg/~matsundf/CorrelationMatrix.m';
filename = 'CorrelationMatrix.m';
if ~exist(filename, 'file')
    fprintf('Downloading %s ...', url);
    websave('CorrelationMatrix.m', url);
    fprintf(' Done!\n');
else
    fprintf('File %s already exists, skipping ...\n', filename);
end