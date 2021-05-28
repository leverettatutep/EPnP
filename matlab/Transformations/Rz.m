function r = Rz(angR,d)
    if nargin == 2
        angR = angR * pi()/180;
    end
    r = eye(3);
    c = cos(angR);
    s = sin(angR);
    r(1,1) = c; r(2,2) = c;
    r(1,2) = -s; r(2,1) = s;
end
