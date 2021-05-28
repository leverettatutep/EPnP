function r = Ry(angR,d)
    if nargin == 2
        angR = angR * pi()/180;
    end
    r = eye(3);
    c = cos(angR);
    s = sin(angR);
    r(1,1) = c; r(3,3) = c;
    r(3,1) = -s; r(1,3) = s;
end
