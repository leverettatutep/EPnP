function [Solution,Alpha,vt,st] = efficient_pnpE(Xw,U,Camera,NumControlPts,...
    Km,hisbetasB,hisscale,hisnumzero,TC2W)

    
    THRESHOLD_REPROJECTION_ERROR=20;%error in degrees of the basis formed by the control points. 

    %define control points in a world coordinate system (centered on the 3d
    %points centroid)
    Cw =define_control_points(NumControlPts);
    CcCor = TC2W * [Cw' ; 1 1 1 1];
    %The ControlPointsInW are arbitrarily chosen

    %compute alphas (linear combination of the control points to represent the 3d
    %points)
    Alpha=compute_alphas(Xw,Cw,NumControlPts);
    numOfPoints = size(Xw,1);

    %Our method
    MaxSolutions = 4;
    Dw = FindDistancesBetweenPoints(Xw);
    [EigVectors,EigVals,vt,st] =findError(Camera,Alpha,U,NumControlPts); 
%     EigVectors = -2 * EigVectors;
    Xc = zeros(numOfPoints,3,MaxSolutions);
    Cc = zeros(NumControlPts,3,MaxSolutions);
    TcTow = zeros(4,4,MaxSolutions);
    Bang = ones(MaxSolutions,1) * 10^15;
    Bdist = ones(MaxSolutions,1) * 10^15;
    error = ones(MaxSolutions,1)*10^15;

    betas = zeros(4,5);
    % Assuming 1 0 eigvalue
    NumZeroEig = 1;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);
    [Bang(NumZeroEig),~,Bdist(NumZeroEig)] = Screws(InvT(TC2W) * TcTow(:,:,NumZeroEig),1);
%     %This also works.
%     %To get the version of Cc as above Cc = cc / scale
%     cc = reshapeVector(EigVectors(:,1));
%     [tctow,scale] = GetScaleTransFromA2B(cc',Cw');

%% This is my version of 2 and 3 zeros    
    %Now assuming 2 0 eigvalue
    NumZeroEig = 2;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);
    [Bang(NumZeroEig),~,Bdist(NumZeroEig)] = Screws(InvT(TC2W) * TcTow(:,:,NumZeroEig),1);

    %Now assuming 3 0 eigvalue
    NumZeroEig = 3;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);
    [Bang(NumZeroEig),~,Bdist(NumZeroEig)] = Screws(InvT(TC2W) * TcTow(:,:,NumZeroEig),1);

    %Now assuming 4 0 eigvalue
    NumZeroEig = 4;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);
    [Bang(NumZeroEig),~,Bdist(NumZeroEig)] = Screws(InvT(TC2W) * TcTow(:,:,NumZeroEig),1);

%% Collecting the solutions
    sol = 1;
    for i=1:MaxSolutions-1
%         if Bang(i) < Bang(sol)
%         if Bdist(i) < Bdist(i)
        if error(i) < error(sol)
            sol = i;
        end
    end
    
    Solution.NumZeros = sol;
    Solution.beta = betas(:,sol);
    Solution.Xc=Xc(:,:,sol);
    Solution.Xw=Xw;
    Solution.Cc=Cc(:,:,sol);
    Solution.Cw=Cw;
    Solution.EigVec=EigVectors;
    Solution.EigVals=EigVals;
    Solution.TC2W=TcTow(:,:,sol);
    Solution.error=error(sol);
    Solution.CcCor = CcCor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)

    %clear all; close all; load reprojection_error_usingRT;
    n=size(Xw,1);

    P=A*[R,T];
    Xw_h=[Xw,ones(n,1)];
    Urep_=(P*Xw_h')';

    %project reference points into the image plane
    Urep=zeros(n,2);
    Urep(:,1)=Urep_(:,1)./Urep_(:,3);
    Urep(:,2)=Urep_(:,2)./Urep_(:,3);

    %reprojection error
    err_=sqrt((U(:,1)-Urep(:,1)).^2+(U(:,2)-Urep(:,2)).^2);
    err=sum(err_)/n;
end

function [Cc, Xc, betas, TcTow, error] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumberOfZeros)
    distC = compute_distancesOfPoints(EigVectors(:,1:NumberOfZeros),Alpha);
    [Cc,Xc,betas]= compute_Cc(EigVectors,distC,Dw,Alpha);
    %todo test
%     [R,T]=getrotT(Xw,Xc);  %solve exterior orientation from C to W
%     TcTow = MakeTFromRTS(R,T);
    [TcTow,aIsBiggerBy,a,b]=GetScaleTransFromA2B(Xc',Xw');
    R = TcTow(1:3,1:3);
    T = TcTow(1:3,4);
    error = reprojection_error_usingRT(Xw,U,R,T,Camera);
end

function CcVectors = GetCcVectors(EigVectors)
    Cc = zeros(4,3,4);
    CcVectors = zeros(3,3,4);
    for i = 1:4
        Cc(:,:,i) = reshape(EigVectors(:,i),4,3);
        for j=1:3
            CcVectors(j,:,i) = Cc(j,:,i) - Cc(4,:,i);
        end
    end
end

function BetaCc = GetBetaCc(CcVectors,Betas)
    for i=1:4
        BetaCc = Betas(i)*CcVectors(:,:,i); 
    end
end
