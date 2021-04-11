% MAIN Illustrates how to use the EPnP algorithm described in:
%
%       Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua.
%       Accurate Non-Iterative O(n) Solution to the PnP Problem. 
%       In Proceedings of ICCV, 2007. 
%
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

clear all; close all;

addpath data;
addpath error;
addpath EPnP;

NumC = 4;
rng(1234);
fprintf('\n---------EPnP--------------\n');
%1.-Generate simulated input data------------------------------------------
for loop = 1:20 %Loop to run multiple times.
load_points=0;
if ~load_points
    %Randomly chooses world points, random camera location, computes image
    %locations of the world points, adds noise to the image values LJE
    n=10; %number of points
    std_noise=10; %noise in the measurements (in pixels)
    std_noise = 0; %LJE
    [Camera,point,tFromCtoW]=generate_noisy_input_data(n,std_noise,'donotplot');
% I plan to define points as 4 by n and 4th row is 1
    TheW = zeros(4,n);
    TheC = TheW;
    for i=1:n
        TheW(:,i) = [point(i).Xworld;1];
        TheC(:,i) = [point(i).Xcam;1];
    end
    TransCtoW = getTransFromA2B(TheW,TheC); %Get T from C to W
    TransWtoC = invT(TransCtoW);
    TheseAre0 = TheC - TransCtoW * TheW;
    save('data\input_data_noise.mat','Camera','point','tFromCtoW','TheW','TheC','TransCtoW','TransWtoC');
else
    load('data\input_data_noise.mat','Camera','point','tFromCtoW','TheW','TheC','TransCtoW','TransWtoC');
    n=size(point,2);
    draw_noisy_input_data(point);
end
% %LJE used to input data to mathematica
% scale = 10000;
% fileID = fopen('uv.csv','w');
% for i=1:50
%     fprintf(fileID,'%d, %d\r\n',round(point(i).Ximg_pix_true*scale,0));
% end
% fclose(fileID);
% %LJE end
%LJE noticed the following are not used
%% 2.-Inputs format--------------------------------
% x3d=zeros(n,4);
% x2d=zeros(n,3); 
Camera=Camera(:,1:3);
%LJE Attaching a 1 at the end of points
for i=1:n
    x3d_h(i,:)=[point(i).Xworld',1]; 
    x2d_h(i,:)=[point(i).Ximg(1:2)',1];

    %world and camera coordinates
    X3d_world(i,:)=point(i).Xworld';
    X3d_cam(i,:)=point(i).Xcam';
end

%% 3.-EPnP----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

%LJE added extra outputs making them available to mathematica
[Bestsol,Mysol,NumOfSV,alphas,Cw,Cc]=efficient_pnp(x3d_h,x2d_h,Camera,NumC);
Rp = Bestsol.R;
Tp = Bestsol.T;
Transp = MakeT(Rp,Tp); %Trans from C to W
Xc = Bestsol.Xc;

Rpe = Mysol.R;
Tpe = Mysol.T;
Transe = MakeT(Rpe,Tpe); %Trans from C to W
Xce = Mysol.Xc;

errorIs = Mysol.error - Bestsol.error;
if abs(errorIs) > 0.1
    display('The number of SV and My error minus his, - means I win  ')
    SvAndError = [NumOfSV errorIs]
end
[Bang,Baxis,Bdist] = Screws(invT(Transp) * Transe,1);
if (abs(Bang)+abs(Bdist)) > .1
    AngleDist = [Bang Bdist]
end
end %of loop
% %LJE writing data to mathematica
% fileID = fopen('alpha.csv','w');
% for i=1:50
%     fprintf(fileID,'%d, %d, %d, %d\r\n',round(scale * alphas(i,:),0));
% end
% fclose(fileID);

% fileID = fopen('cw.csv','w');
% for i=1:4
%     fprintf(fileID,'%d, %d, %d\r\n',round(scale * Cw(i,:),0));
% end
% fclose(fileID);
% 
% fileID = fopen('cc.csv','w');
% for i=1:4
%     fprintf(fileID,'%d, %d, %d\r\n',round(scale * Cc(i,:),0));
% end
% fclose(fileID);

%draw Results
for i=1:n
    point(i).Xcam_est=Xc(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP (Old)',h);
xlim([-2 2]); ylim([-2 2]);

%compute error
error=reprojection_error_usingRT(Xw,U,Rp,Tp,Camera);
fprintf('error EPnP: %.3f\n',error);


%3.-EPnP_GAUSS_NEWTON----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

[Rp,Tp,Xc,sol]=efficient_pnp_gauss(x3d_h,x2d_h,Camera);

%draw Results
for i=1:n
    point(i).Xcam_est=Xc(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP Gauss Newton',h);

%compute error
error=reprojection_error_usingRT(Xw,U,Rp,Tp,Camera);
fprintf('error EPnP_Gauss_Newton: %.3f\n',error);
xlim([-2 2]); ylim([-2 2]);




