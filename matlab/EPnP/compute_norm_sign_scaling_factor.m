function [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw)
 
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


n=size(Xw,1); %number of data points

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
numRows = size(X1,1);
if numRows == 12
    Cc_ = reshape(X1,[3,4])'; %Cc row 1 is CP 1, Cc row 2 is CP 2 
else
    Cc_ = reshape(X1,[3,3])';
end

%position of reference points in camera coordinates
Xc_=Alph*Cc_;

%compute distances in world coordinates w.r.t. the centroid
dist_w=FindDistancesBetweenPoints(Xw);

%compute distances in camera coordinates w.r.t. the centroid
dist_c=FindDistancesBetweenPoints(Xc_);
 
%least squares solution for the scale factor
%If there was a single point then dist_w and dist_c would be 1 by 1
%Solution would be dist_c beta = dist_W and beta = inv(dist_c) * dist_w
%So for multiple points (n) there are n equations like those above
%There are n-1 too many equations so you must find the BEST or least square
%Solution which is: A x = B -> A' * A x = A' B -> x = inv(A' * A)(A' * B)
sc=1/(inv(dist_c'*dist_c)*dist_c'*dist_w); %scale is 1/beta

%scale position of the control points
Cc=Cc_/sc;

%rescaled position of the reference points
Xc=Alph*Cc;

%change the sign if necessary. z negative is no possible in camera
%coordinates
neg_z=find(Xc(:,3)<0);
if size(neg_z,1)>=1
    sc=-sc;
    Xc=Xc*(-1);
end
