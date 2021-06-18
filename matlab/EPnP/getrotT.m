function [R, T]=getrotT(wpts,cpts);
  
% This routine solves the exterior orientation problem for a point cloud
%  given in both camera and world coordinates. 
  
% wpts = 3D points in arbitrary reference frame
% cpts = 3D points in camera reference frame
  
    n=size(wpts,1);
    M=zeros(3);

    ccent=mean(cpts);
    wcent=mean(wpts);

%This changes the origins so they are measured from the centers
%this makes sure these points are off by a rotation and not a 
%translation
    for i=1:3
      cpts(:,i)=cpts(:,i)-ccent(i)*ones(n,1);
      wpts(:,i)=wpts(:,i)-wcent(i)*ones(n,1);
    end

%Not sure why it is written this way this is just
%M = cpts' * wpts
    for i=1:n
       M=M+cpts(i,:)'*wpts(i,:);
    end

    [U S V]=svd(M);
    R=U*V';
%There is really no reason why the R matrix must be right handed
%Left handed frames have negative determinates. So this ensures
%it is right handed.
    if det(R)<0
      R=-R;
    end

%This allows for the translation between the point centers.
%I believe I have work that shows this will be the optimal translation.
    T=ccent'-R*wcent';
end

