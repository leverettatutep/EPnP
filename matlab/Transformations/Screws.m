function [ang, axis, dist, scale] = Screws(t,ASecondInputMeansReturnDegrees)
    scale = nthroot(det(t),3);
    t = t/scale;
    axis = zeros(3,1);
    cosine = (t(1,1)+t(2,2)+t(3,3)-1)/2;
    kysine = t(1,3)-t(3,1);
    kzsine = t(2,1)-t(1,2);
    kxsine = t(3,2)-t(2,3);
    sine = sqrt((kysine^2+kzsine^2+kxsine^2)/4);
    if nargin > 1
        ang = atan2d(sine,cosine);
    else
        ang = atan2(sine,cosine);
    end
    if sine == 0
        axis(1)=1;
    else
        axis =[kxsine/(2*sine);
               kysine/(2*sine);
               kzsine/(2*sine)];
    end
    dist = 0;
    if size(t,1) == 4
        dist = sqrt(t(1:3,4)' * t(1:3,4));
    end
%     dist = dist / scale;
end
    