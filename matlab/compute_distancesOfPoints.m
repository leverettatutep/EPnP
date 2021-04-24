function dist_c=compute_distancesOfPoints(EigenVector,Alph)
    NumToCompute = size(EigenVector,2);
    dist_c = zeros(size(Alph,1),NumToCompute);
    for i=1:NumToCompute
        Cc_ = reshapeVector(EigenVector(:,i));
        Xc_=Alph*Cc_;
        dist_c(:,i) = FindDistancesBetweenPoints(Xc_);
    end
end
