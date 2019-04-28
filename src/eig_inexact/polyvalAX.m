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

function Y = polyvalAX(p,A,X,a,b)                                               
% Y = polyvalAX(p,A,X,a,b)                                                      
                                                                                
degree = numel(p) - 1;                                                          
if degree == 0; Y = p*X; return; end                                            
                                                                                
c0 = (a + b)/(a - b);                                                           
c1 = 2 / (b - a);                                                               
B = @(X) c0*X + c1*A(X);                                                        
                                                                                
Y = p(1)*X; del = eps;                                                          
for j = 1:degree                                                                
    %Y = B(Y) + p(j+1)*X                                                        
    Y = B(Y);                                                                   
    if abs(p(j+1)) > del                                                        
        Y = Y + p(j+1)*X;                                                       
    end                                                                         
end                                                                             
end
