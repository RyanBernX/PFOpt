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

%% uninstall_low_rank_opt.m
%  Remove all low-rank-opt files

if exist('low-rank-opt', 'dir')
    rmdir('low-rank-opt', 's');
end

if exist('low-rank-opt.zip', 'file')
    delete('low-rank-opt.zip');
end