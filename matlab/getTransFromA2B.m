function trans = getTransFromA2B(a,b)
%a are points 4 by n
%b are points 4 by n
    %His program expects points are n by 3
    [R,T]=getrotT(a(1:3,:)',b(1:3,:)'); %Get T from W to C
    trans = eye(4);
    trans(1:3,1:3) = R;
    trans(1:3,4) = T;
end
