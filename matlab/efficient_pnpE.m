function [Solution,Alpha,vt,st] = efficient_pnpE(Xw,U,Camera,NumControlPts,...
    Km,hisbetasB,hisscale,hisnumzero)
    THRESHOLD_REPROJECTION_ERROR=20;%error in degrees of the basis formed by the control points. 

    %define control points in a world coordinate system (centered on the 3d
    %points centroid)
    Cw =define_control_points(NumControlPts);
    %The ControlPointsInW are arbitrarily chosen

    %compute alphas (linear combination of the control points to represent the 3d
    %points)
    Alpha=compute_alphas(Xw,Cw,NumControlPts);
    numOfPoints = size(Xw,1);

    %Our method
    MaxSolutions = 5; %4;
    Dw = FindDistancesBetweenPoints(Xw);
    [EigVectors,EigVals,vt,st] =findError(Camera,Alpha,U,NumControlPts); 
%     EigVectors = -2 * EigVectors;
    Xc = zeros(numOfPoints,3,MaxSolutions);
    Cc = zeros(NumControlPts,3,MaxSolutions);
    TcTow = zeros(4,4,MaxSolutions);
    error = ones(MaxSolutions,1)*10^15;

    betas = zeros(4,5);
    % Assuming 1 0 eigvalue
    NumZeroEig = 1;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);

%% This is my version of 2 and 3 zeros    
    %Now assuming 2 0 eigvalue
    NumZeroEig = 2;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);

    %Now assuming 3 0 eigvalue
    NumZeroEig = 3;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);

    %Now assuming 4 0 eigvalue
    NumZeroEig = 4;
    [Cc(:,:,NumZeroEig), Xc(:,:,NumZeroEig), ...
        betas(:,NumZeroEig), TcTow(:,:,NumZeroEig), ...
        error(NumZeroEig)] = ComputeSol(...
        EigVectors, Alpha, Dw, Xw, U, Camera, NumZeroEig);

    %% This is my emulation of his betas.
    hisbetas = zeros(1,4);
    if hisnumzero == 2
        hisbetas(2) = hisbetasB(1);
        hisbetas(1) = hisbetasB(2);
        hisbetas = hisbetas/hisscale;
    end
    if hisnumzero == 3
        hisbetas(3) = hisbetasB(1);
        hisbetas(2) = hisbetasB(2);
        hisbetas(1) = hisbetasB(3);
        hisbetas = hisbetas/hisscale;
    end
    if hisnumzero == 4
        hisbetas(4) = hisbetasB(1);
        hisbetas(3) = hisbetasB(2);
        hisbetas(2) = hisbetasB(3);
        hisbetas(1) = hisbetasB(4);
        hisbetas = hisbetas/hisscale;
    end
    betas(:,5) = hisbetas;
    if hisnumzero > 1
        for i=1:4 %make sure they are correct sign
            if EigVectors(1,i) * Km(1,i) < 0
                EigVectors(:,i) = -EigVectors(:,i);
            end
        end
    end
    
    %Now use his betas with my eigenvector
    distChis = compute_distancesOfPoints(EigVectors(:,1:4),Alpha);
    [Cc(:,:,5),Xc(:,:,5),~]= ...
        compute_Cc(EigVectors,distChis,Dw,Alpha, betas(:,5));
    [R,T]=getrotT(Xw,Xc(:,:,5));  %solve exterior orientation from C to W
    TcTow(:,:,5) = MakeT(R,T);
    error(5) =reprojection_error_usingRT(Xw,U,R,T,Camera);

%% Collecting the solutions
    sol = 1;
    for i=1:MaxSolutions-1
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R, T]=getrotT(wpts,cpts)

    % This routine solves the exterior orientation problem for a point cloud
    %  given in both camera and world coordinates. 

    % wpts = 3D points in arbitrary reference frame
    % cpts = 3D points in camera reference frame

    n=size(wpts,1);
    M=zeros(3);

    ccent=mean(cpts);
    wcent=mean(wpts);

    for i=1:3
      cpts(:,i)=cpts(:,i)-ccent(i)*ones(n,1);
      wpts(:,i)=wpts(:,i)-wcent(i)*ones(n,1);
    end
    for i=1:n
       M=M+cpts(i,:)'*wpts(i,:);
    end
    [U S V]=svd(M);
    R=U*V';
    if det(R)<0
      R=-R;
    end
    T=ccent'-R*wcent';
    % 
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
    [R,T]=getrotT(Xw,Xc);  %solve exterior orientation from C to W
    TcTow = MakeT(R,T);
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

function error = NormalVectorError(CcVectors)