function scale =FindDistances(Cc,Cw)

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
centr_Cc = mean(Cc);
centroid_Cc=repmat(centr_Cc,[4,1]);
tmp1=Cc-centroid_Cc;
dist_Cc=sqrt(sum(tmp1.^2,2));

centr_Cw = mean(Cw);
centroid_Cw=repmat(centr_Cw,[4,1]);
tmp1=Cw-centroid_Cw;
dist_Cw=sqrt(sum(tmp1.^2,2));


% %compute distances in world coordinates w.r.t. the centroid
% centr_w=mean(Xw);
% centroid_w=repmat(centr_w,[n,1]);
% tmp1=Xw-centroid_w;
% dist_w=sqrt(sum(tmp1.^2,2));

% %compute distances in camera coordinates w.r.t. the centroid
% centr_c=mean(Xc_);
% centroid_c=repmat(centr_c,[n,1]);
% tmp2=Xc_-centroid_c;
% dist_c=sqrt(sum(tmp2.^2,2));
 
%least squares solution for the scale factor
%Let Noah deal with this
sc=1/(inv(dist_Cc'*dist_Cc)*dist_Cc'*dist_Cw);

%scale position of the control points
Cc=Cc/sc;
scale = 1/sc;
% %rescaled position of the reference points
% Xc=Alph*Cc;
% 
% %change the sign if necessary. z negative is no possible in camera
% %coordinates
% neg_z=find(Xc(:,3)<0);
% if size(neg_z,1)>=1
%     sc=-sc;
%     Xc=Xc*(-1);
% end

end


        
        
