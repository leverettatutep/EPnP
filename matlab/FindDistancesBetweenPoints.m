function distance = FindDistancesBetweenPoints(pts)
    n=size(pts,1);    %pts is n by 3 and n = number of points
    centr_w=mean(pts); %Find the average point (the center), 1 by 3
    centroid_w=repmat(centr_w,[n,1]); %Copy the center into n rows by 3
    tmp1=pts-centroid_w; %n points from center to pts
    distance=sqrt(sum(tmp1.^2,2)); %Each term in the n by 3 is squared 
            %then columns are summed forming an n by 1 vector which are
            %lengths of each point
end