function [theVs, theVals, vt, st] = findError(camera, alpha, uv, NumC)
    n = size(alpha,1);
    uv = uv';
    NumUnk = NumC * 3;
    error = zeros(2,NumUnk);
    equation = zeros(NumUnk,NumUnk);
    for i =1:n
        errorI = GetErrorI(camera,alpha(i,:),uv(:,i),NumC);
        equation = equation + GetEquationI(camera,alpha(i,:),uv(:,i),errorI);
        error = error + errorI;
    end
    %A x = B
    A = equation;
    [V,S] = eig(A);
    [vt,st] = eig(A * A');
    st = sqrt(st); %These are the same as S
    %vt is +or- V
%The 12 unknowns are listed as [Cx1 Cx2 Cx3 Cx4 Cy1 Cy2 Cy3 Cy4 Cz1
     %Cz2 Cz3 Cz4]
%However the original paper has them as:
%Cx1 Cy1 Cz1; Cx2 Cy2 Cz2 ...
    V1 = reorderV(V(:,1)',NumC);
    V2 = reorderV(V(:,2)',NumC);
    V3 = reorderV(V(:,3)',NumC);
    V4 = reorderV(V(:,4)',NumC);
    theVs = [V1 V2 V3 V4];
    theVals = [S(1,1) S(2,2) S(3,3) S(4,4)];
%      if NumC == 3
%          Cc = reshape(V1,[3,3]);
%          V1 = reshape(Cc',[9,1]);
%          Ccp = [Cc' ; 1 1 1];
%      else
%          Cc = reshape(V1,[4,3]);
%          V1 = reshape(Cc',[12,1]);
%          Ccp = [Cc' ; 1 1 1 1];
%      end

end

function vector = GetVector(c,tail, head)
    vector = c(:,head)-c(:,tail);
end
function leng = GetLength(vector)
    sqrt(vector' * vector);
end

% {fu (xc[1] \[Alpha][1] + xc[2] \[Alpha][2] + xc[3] \[Alpha][3] + 
%      xc[4] \[Alpha][4]) + (uc - ui) (zc[1] \[Alpha][1] + 
%      zc[2] \[Alpha][2] + zc[3] \[Alpha][3] + zc[4] \[Alpha][4]), 

%  fv (yc[1] \[Alpha][1] + yc[2] \[Alpha][2] + yc[3] \[Alpha][3] + 
%      yc[4] \[Alpha][4]) + (vc - vi) (zc[1] \[Alpha][1] + 
%      zc[2] \[Alpha][2] + zc[3] \[Alpha][3] + zc[4] \[Alpha][4])}
function equationI = GetEquationI(fuv2,alphaI,uvI,errorI)
     %The 12 unknowns are listed as [Cx1 Cx2 Cx3 Cx4 Cy1 Cy2 Cy3 Cy4 Cz1
     %Cz2 Cz3 Cz4]
     fu = fuv2(1,1);
     fv = fuv2(2,2);
     uc = fuv2(1,3);
     vc = fuv2(2,3);
%      fuAlphaI = fu * alphaI';
%      fvAlphaI = fv * alphaI';
     ucmuI = (uc - uvI(1)) * alphaI;
     vcmvI = (vc - uvI(2)) * alphaI;
%      for i=1:4
         equationI = [2 * fu * alphaI' * errorI(1,:); 
             2 * fv * alphaI' * errorI(2,:);
             2 * ucmuI' * errorI(1,:) + ...
             2 * vcmvI' * errorI(2,:) ];
%      end
end

function errorI = GetErrorI(fuv2,alphaI,uvI,NumC)
     %The 12 unknowns are listed as [Cx1 Cx2 Cx3 Cx4 Cy1 Cy2 Cy3 Cy4 Cz1
     %Cz2 Cz3 Cz4]
     fu = fuv2(1,1);
     fv = fuv2(2,2);
     uc = fuv2(1,3);
     vc = fuv2(2,3);
     fuAlphaI = fu * alphaI;
     fvAlphaI = fv * alphaI;
     ucmuI = (uc - uvI(1)) * alphaI;
     vcmvI = (vc - uvI(2)) * alphaI;
     if NumC == 3
        errorU = [fuAlphaI 0 0 0 ucmuI];
        errorV = [0 0 0 fvAlphaI vcmvI];
     else
        errorU = [fuAlphaI 0 0 0 0 ucmuI];
        errorV = [0 0 0 0 fvAlphaI vcmvI];
     end         
     errorI = [errorU; errorV];
end
     