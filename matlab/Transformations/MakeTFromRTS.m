function Trans = MakeT(r,t,s)
    if nargin < 3
        s = 1;
    end
    Trans = eye(4);
    Trans(1:3,1:3)= r(1:3,1:3);
    Trans(1:3,4) = t(1:3);
    Trans = s * Trans;
    Trans(4,4) = 1;
end
