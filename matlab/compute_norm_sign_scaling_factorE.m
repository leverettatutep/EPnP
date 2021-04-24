function [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw)
 
    numRows = size(X1,1);
    if numRows == 12
        Cc_ = reshape(X1,[3,4])'; %Cc row 1 is CP 1, Cc row 2 is CP 2 
    else
        Cc_ = reshape(X1,[3,3])';
    end

    Xc_=Alph*Cc_;

    %compute distances in world coordinates w.r.t. the centroid
    dist_c=FindDistancesBetweenPoints(Xc_);

end
