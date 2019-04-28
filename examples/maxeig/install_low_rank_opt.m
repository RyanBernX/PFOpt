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

%% install_low_rank_opt.m
%  Download the most recent low-rank-opt from
%  https://www.cs.ubc.ca/~mpf/pubs/low-rank-spectral-optimization-via-gauge-duality/
%  and setup the examples

filename = 'low-rank-opt.zip';
url = 'http://www.cs.ubc.ca/~mpf/low-rank-opt/low-rank-opt.zip';

if ~exist(filename, 'file')
    fprintf('Downloading %s ...', url);
    websave('low-rank-opt.zip', url);
    fprintf(' Done!\n');
else
    fprintf('File %s already exists, skipping ...\n', filename);
end

if ~exist('low-rank-opt', 'dir')
    fprintf('Unzipping %s ...', filename);
    unzip('low-rank-opt.zip');
    fprintf(' Done!\n');
else
    fprintf('Folder low-rank-opt already exists, skipping ...\n');
end

fprintf('Installing PFGAUGE ...');
copyfile('low-rank-opt-modified', 'low-rank-opt');
fprintf(' Done!\n');

fprintf('\nPFGAUGE has been installed successfully.\n');
