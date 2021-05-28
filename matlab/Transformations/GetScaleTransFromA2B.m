function [trans,aIsBiggerBy,a,b] = GetScaleTransFromA2B(a,b)
% points a may have a different scale than b. 
% First it determines an optimal R to move the a frame to the b.
% If the R is left handed (det < 0) you take - R this is essentially the
% same R but the unit vectors flip so it becomes right handed. Then it
% computes an optimal scale. Then it determines a displacement. The return
% trans is an ordinary transformation, it does NOT have the scale in it.
% The return a is the input a scaled so it matches the b scale. The output
% b is the input b scaled so it matches the a scale. To use the trans
% correctly you should do this: PointsInA = trans * PointsInB and the
% points in A will have the same scale as the scale of B. 
% The points must be in columns. 
%Ensure they have 3 rows, if not there is an error
    n=size(a,2);
    ar = size(a,1);
    br = size(b,1);
    if (ar < 3) || (ar >4) || (br < 3) || (br > 4)
        error('The inputs must be points in columns so rows must be at least 3');
    end
    Bpts = b(1:3,:)';
    Apts = a(1:3,:)';

    acent=mean(Apts,1);
    bcent=mean(Bpts,1);

    for i=1:3
      apts(:,i)=Apts(:,i)-acent(i)*ones(n,1);
      bpts(:,i)=Bpts(:,i)-bcent(i)*ones(n,1);
    end
    
    M=zeros(3);
    for i=1:n
       M=M+apts(i,:)'*bpts(i,:);
    end

    [U, S, V]=svd(M);
    R=U*V';
    if det(R)<0
      R=-R;
    end
    
    scale = GetNumer(bpts,apts,R)/GetDenom(apts);
    T= -((bcent')-R*acent'*scale);
    T = acent'*scale - R*bcent';
%     if scale < 0
%         scale = - scale;
%     end
    trans = eye(4);
    trans(1:3,1:3) = R;
    trans(1:3,4) = T;
%     trans = trans * scale;
    trans(4,4)=1;
    aIsBiggerBy = 1/scale;
    a = trans(1:3,1:3)*Bpts'+trans(1:3,4);
    b = trans(1:3,1:3)'*(Apts'-trans(1:3,4));
end

function denom = GetDenom(y)
    n = size(y,1);
    denom = 0;
    for i=1:n
        denom = denom + y(i,:) * y(i,:)';
    end
end

function numer = GetNumer(x,y,r)
    n = size(y,1);
    numer = 0;
    for i=1:n
        numer = numer + x(i,:) * r' * y(i,:)';
    end
end
