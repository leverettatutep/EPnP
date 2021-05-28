function [Cc,Xc,betas]=compute_Cc(Vectors,Dc,Dw,Alpha,hisB)

%least squares solution for the scale factor
%If there was a single point then dist_w and dist_c would be 1 by 1
%Solution would be dist_c beta = dist_W and beta = inv(dist_c) * dist_w
%So for multiple points (n) there are n equations like those above
%There are n-1 too many equations so you must find the BEST or least square
%Solution which is: A x = B -> A' * A x = A' B -> x = inv(A' * A)(A' * B)
betas = zeros(4,1);
if nargin < 5
    beta=(inv(Dc'*Dc)*Dc'*Dw); 
else
    beta = hisB;
end
CcVector = zeros(size(Vectors,1),1);
if size(beta,1) == 1
    CcVector = beta * Vectors(:,1);
    betas(1,1) = beta;
else
    for i=1:size(beta,1)
        CcVector = CcVector(:) + beta(i)*Vectors(:,i);
        betas(i,1) = beta(i);
    end
end
Cc = reshapeVector(CcVector);

%rescaled position of the reference points
Xc=Alpha*Cc;

%change the sign if necessary. z negative is no possible in camera
%coordinates
neg_z=find(Xc(:,3)<0);
if size(neg_z,1)>=1
    Xc=Xc*(-1);
    Cc = Cc * (-1);
end
