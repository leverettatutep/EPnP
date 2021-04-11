function Trans = MakeT(r,t)
    Trans = eye(4);
    Trans(1:3,1:3)= r(1:3,1:3);
    Trans(1:3,4) = t(1:3);
end
