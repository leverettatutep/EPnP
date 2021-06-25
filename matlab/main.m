clear all; close all; clc;

addpath data;
addpath error;
addpath EPnP;
addpath Transformations;

%Parameters to control
NumC = 4;
NumPoints = 10;
Noise = 10;
% Noise = 0;
Loops = 1000;
load_points=0;

winners = zeros(Loops,1);
% numOfZeros = zeros(Loops,1);
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
    [Camera,point,tFromCtoW,centroid]=generate_noisy_input_data(n,std_noise,'donotplot');
% I plan to define points as 4 by n and 4th row is 1
    TheW = zeros(4,n);
    TheC = TheW;
    for i=1:n
        TheW(:,i) = [point(i).Xworld;1];
        TheC(:,i) = [point(i).Xcam;1];
    end
    [TransCtoW,CIsBiggerBy] = GetScaleTransFromA2B(TheC,TheW); %Get T from C to W
    [TransWtoC,WIsBiggerBy] = InvT(TransCtoW);
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
TransHis = MakeTFromRTS(Rp,Tp); %Trans from C to W same as Bestsol.TC2W
Xc = Bestsol.Xc;
pCw = MakeHomoPts(Cw');
pCc = MakeHomoPts(-Cc');

% max(max(Bestsol.TC2W * MakeHomoPts(Bestsol.Xw') - MakeHomoPts(Bestsol.Xc')))
% max(max(Bestsol.TC2W * pCw - pCc))
[TC2Whis,scaleHis]=GetScaleTransFromA2B(pCc,pCw);
% max(max(TC2Whis * inv(GetScaleTransFromA2B(-Cc',Cw')) - eye(4)))
% max(max(GetScaleTransFromA2B(Bestsol.Xc',Bestsol.Xw')*inv(TC2Whis) - eye(4)))
% TC2Whis

[Bang,~,Bdist] = Screws(InvT(tFromCtoW) * TransHis,1);
TransErrorh = sqrt(Bang^2+Bdist^2);
him(loop,:) = [Bestsol.error Bestsol.NumZeros Bang Bdist];

%% My solution-------------------------------------------------
%LJE added extra outputs making them available to mathematica
Weights = ones(NumPoints,1);
[Mysol,alphas,vt,st]=efficient_pnpE(Xw,U,Camera,NumC,Weights,Bestsol.EigVec,Bestsol.beta,Bestsol.scale,Bestsol.NumZeros,TransCtoW);
Xce = Mysol.Xc;
Xwe = Mysol.Xw;
Cce = Mysol.Cc;
Cwe = Mysol.Cw;
EigVec = Mysol.EigVec;
EigVals = Mysol.EigVals/2; %My vals are two times his
TC2W = Mysol.TC2W;
error = Mysol.error;
% if Mysol.NumZeros > 1
%     Mysol.NumZeros
%     Cce
% end

% max(max(Bestsol.TC2W * inv(Mysol.TC2W) - eye(4)))
% max(max(Bestsol.Xc - Mysol.Xc))
pCwe = MakeHomoPts(Mysol.Cw');
pCce = MakeHomoPts(-Mysol.Cc');

% max(max(Mysol.TC2W * MakeHomoPts(Mysol.Xw') - MakeHomoPts(Mysol.Xc')))
% max(max(Mysol.TC2W * pCwe - pCce))
[TC2We,scalee]=GetScaleTransFromA2B(pCce,pCwe);
% max(max(TC2We * inv(GetScaleTransFromA2B(-Cce',Cwe')) - eye(4)))
% max(max(GetScaleTransFromA2B(Mysol.Xc',Mysol.Xw')*inv(TC2We) - eye(4)))
% max(max(TC2We * inv(TC2W) - eye(4)))

[Bang,~,Bdist] = Screws(InvT(tFromCtoW) * TC2W,1);
TransError = sqrt(Bang^2+Bdist^2);
me(loop,:) = [error Mysol.NumZeros Bang Bdist];

%% Comparing ----------------------------------------------------
% errorIs = error - Bestsol.error;
% eigDif(loop,:) = [max(abs(EigVals - Bestsol.EigVals')) max(max(abs(abs(EigVec) - abs(Bestsol.EigVec))))];
%my eigenvectors are +or- his
% The eigenvalues are identical always
% The eigenvectors can differ
    loop
end %of loop

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

HimMe = him - me;




