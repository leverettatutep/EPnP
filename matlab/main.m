clear all; close all; clc;

addpath data;
addpath error;
addpath EPnP;

%Parameters to control
NumC = 4;
NumPoints = 10;
Noise = 10;
Noise = 10;
Loops = 1000;
load_points=0;

winners = zeros(Loops,1);
numOfZeros = zeros(Loops,1);
me = zeros(Loops,4);
him = zeros(Loops,4);
eigDif = zeros(Loops,2);

rng(1234);
fprintf('\n---------EPnP--------------\n');
%%1.-Generate simulated input data------------------------------------------
for loop = 1:Loops %Loop to run multiple times.
if ~load_points
    %Randomly chooses world points, random camera location, computes image
    %locations of the world points, adds noise to the image values LJE
    n=NumPoints; %number of points
    std_noise=Noise; %noise in the measurements (in pixels)
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
    TheseAre0 = max(max(abs(TheC - TransCtoW * TheW)));
    [Bang,~,Bdist] = Screws(tFromCtoW * TransWtoC,1);
    TheseAreAlso0 = [abs(Bang),abs(Bdist)];
    if abs(Bang)+abs(Bdist)+TheseAre0 > 0.0001
        display('There is an error.');
    end
%     if (abs(Bang)+abs(Bdist)) > .1
%         InitialAngleDist = [Bang Bdist]
%     end
%      save('data\input_data_noise.mat','Camera','point','tFromCtoW','TheW','TheC','TransCtoW','TransWtoC');
     fclose('all');
else
    load('data\input_data_noise.mat','Camera','point','tFromCtoW','TheW','TheC','TransCtoW','TransWtoC');
    n=size(point,2);
    draw_noisy_input_data(point);
end
%% 2.-Inputs format--------------------------------
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

%% His solution-------------------------------------------------
[Bestsol,alphas,Cw,Cc]=efficient_pnp(x3d_h,x2d_h,Camera,NumC);
Rp = Bestsol.R;
Tp = Bestsol.T;
TransHis = MakeT(Rp,Tp); %Trans from C to W
Xc = Bestsol.Xc;

[Bang,~,Bdist] = Screws(invT(tFromCtoW) * TransHis);
TransErrorh = sqrt(Bang^2+Bdist^2);
him(loop,:) = [Bestsol.error Bestsol.NumZeros Bang Bdist];

%% My solution-------------------------------------------------
%LJE added extra outputs making them available to mathematica
[Mysol,alphas,vt,st]=efficient_pnpE(Xw,U,Camera,NumC,Bestsol.EigVec,Bestsol.beta,Bestsol.scale,Bestsol.NumZeros);
Xce = Mysol.Xc;
Xwe = Mysol.Xw;
Cce = Mysol.Cc;
Cwe = Mysol.Cw;
EigVec = Mysol.EigVec;
EigVals = Mysol.EigVals/2; %My vals are two times his
TC2W = Mysol.TC2W;
error = Mysol.error;
if Mysol.NumZeros > 1
    Mysol.NumZeros
    Cce
end

[Bang,~,Bdist] = Screws(invT(tFromCtoW) * TC2W);
TransError = sqrt(Bang^2+Bdist^2);
me(loop,:) = [error Mysol.NumZeros Bang Bdist];

%% Comparing ----------------------------------------------------
errorIs = error - Bestsol.error;
if errorIs < 38 && errorIs > 24.3
    display ('here now');
end
if (TransErrorh - TransError) < -1.88
    display ('this is it');
end
eigDif(loop,:) = [max(abs(EigVals - Bestsol.EigVals')) max(max(abs(abs(EigVec) - abs(Bestsol.EigVec))))];
%my eigenvectors are +or- his
% The eigenvalues are identical always
% The eigenvectors can differ
if error > Bestsol.error
    display('He wins')
else
    display('I win')
end

end %of loop
% [max(winners) min(winners) max(numOfZeros)]
% sumMe = 0;
% sumHe = 0;
% numMe = 0;
% numHe = 0;
% for loop=1:Loops
%     if winners(loop) < 0
%         numMe = numMe+1;
%         sumMe = sumMe + winners(loop);
%     else
%         numHe = numHe+1;
%         sumHe = sumHe + winners(loop);
%     end
% end
% [sumMe/numMe sumHe/numHe]

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


%%3.-EPnP_GAUSS_NEWTON----------------------------------------------------
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




