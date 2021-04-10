function distance = FindDistancesBetweenPoints(pts)
    n=size(pts,1);
    centr_w=mean(pts);
    centroid_w=repmat(centr_w,[n,1]);
    tmp1=pts-centroid_w;
    distance=sqrt(sum(tmp1.^2,2));
end