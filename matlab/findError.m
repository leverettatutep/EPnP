function Cc = findError(camera, alpha, uv)
    n = size(alpha,1);
    uv = uv';
%     c = [1 1 1; 2 1 1; 1 2 1; 1 1 2]';
%     cans = [c(1,:) c(2,:) c(3,:)]';
%     camera = [800 0 320; 0 800 240; 0 0 1];
%     alpha = rand(n,4);
%     wuv = camera * c * (alpha');
%     uv = zeros(2,n);
%     wi = wuv(3,:);
%     for i=1:n
%         uv(:,i) = wuv(1:2,i)/wi(i);
%     end

    error = zeros(2,12);
    equation = zeros(12,12);
    for i =1:n
        errorI = GetErrorI(camera,alpha(i,:),uv(:,i));
        equation = equation + GetEquationI(camera,alpha(i,:),uv(:,i),errorI);
        error = error + errorI;
    end
%     equation * cans;
%     error * cans;
    %A x = B
    A = equation;
    [V,S] = eig(A);
    %Expect only one zero eigenvalue
    V1 = V(:,1)';
%The 12 unknowns are listed as [Cx1 Cx2 Cx3 Cx4 Cy1 Cy2 Cy3 Cy4 Cz1
     %Cz2 Cz3 Cz4]
     Cc = [V1(1:4);V1(5:8);V1(9:12)]';
     
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

function errorI = GetErrorI(fuv2,alphaI,uvI)
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
     errorU = [fuAlphaI 0 0 0 0 ucmuI];
     errorV = [0 0 0 0 fvAlphaI vcmvI];
     errorI = [errorU; errorV];
end
     