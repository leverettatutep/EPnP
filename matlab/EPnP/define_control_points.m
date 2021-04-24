function Cw=define_control_points(number)

% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

if nargin<1 
    number=4; 
end
if number == 3
    Cw = eye(3);
else
    Cw=[1 0 0;
        0 1 0;
        0 0 1;
        0 0 0];
end
Cw = 1 * Cw;
%LJE testing to see impact of other values of Cw
% Cw=[1 0 0;
%     0 1 0;
%     0 0 1;
%     1 1 1];
