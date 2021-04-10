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
%     std_noise = 0; %LJE
    [A,point,Rt,centroid]=generate_noisy_input_data(n,std_noise,'donotplot');
    Trans = Rt;
    Trans(1:4,4) = -Trans * [centroid ; 1];
    Trans(4,4) = 1;
    this = Trans * [(point(1).Xcam);1];
    this - ([point(1).Xworld;1]);
    save('data\input_data_noise.mat','A','point','Rt','Trans');
else
    load('data\input_data_noise.mat','A','point','Rt','Trans');
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
A=A(:,1:3);
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
[Bestsol,Mysol,NumOfSV,alphas,Cw,Cc]=efficient_pnp(x3d_h,x2d_h,A,Trans,NumC);
Rp = Bestsol.R;
Tp = Bestsol.T;
Xc = Bestsol.Xc;
errorIs = Mysol.error - Bestsol.error;
if abs(errorIs) > 0.1
    display('The number of SV and My error minus his, - means I win  ')
    NumOfSV
    errorIs
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
error=reprojection_error_usingRT(Xw,U,Rp,Tp,A);
fprintf('error EPnP: %.3f\n',error);


%3.-EPnP_GAUSS_NEWTON----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

[Rp,Tp,Xc,sol]=efficient_pnp_gauss(x3d_h,x2d_h,A);

%draw Results
for i=1:n
    point(i).Xcam_est=Xc(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP Gauss Newton',h);

%compute error
error=reprojection_error_usingRT(Xw,U,Rp,Tp,A);
fprintf('error EPnP_Gauss_Newton: %.3f\n',error);
xlim([-2 2]); ylim([-2 2]);




