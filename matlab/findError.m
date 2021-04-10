function V1 = findError(camera, alpha, uv, Trans, NumC)
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
%These lines are intended to investigate the number of singular values.
    sizes = zeros(1,NumUnk);
    for col=1:NumUnk
        sizes(col) = norm(S(:,col));%/norm(S(:,NumUnk));
    end
    sizes;
%end of the SV section
    %Expect only one zero eigenvalue
    V1 = V(:,1)';
%The 12 unknowns are listed as [Cx1 Cx2 Cx3 Cx4 Cy1 Cy2 Cy3 Cy4 Cz1
     %Cz2 Cz3 Cz4]
%However the original paper has them as:
%Cx1 Cy1 Cz1; Cx2 Cy2 Cz2 ...
     if NumC == 3
         Cc = reshape(V1,[3,3]);
         V1 = reshape(Cc',[9,1]);
         Ccp = [Cc' ; 1 1 1];
     else
         Cc = reshape(V1,[4,3]);
         V1 = reshape(Cc',[12,1]);
         Ccp = [Cc' ; 1 1 1 1];
     end
     Trans * Ccp;
%     S
%     V
%     AtA=A'*A;
%     [V,S]=eig(AtA);
%     S
%     V
% 
%     z = 11;
%     A8 = A(1:z,1:z);
%     [V,S,W] = eig(A8);
%     S
%     V
% 
%     AtA8=A8'*A8;
%     [V,S]=eig(AtA8);
%     S
%     V

%     B = zeros(12,1);
%     B8 = B(1:z);
%     A8\B8
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
     