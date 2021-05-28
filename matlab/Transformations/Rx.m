function r = Rx(angR,d)
    if nargin == 2
        angR = angR * pi()/180;
    end
    r = eye(3);
    c = cos(angR);
    s = sin(angR);
    r(2,2) = c; r(3,3) = c;
    r(3,2) = s; r(2,3) = -s;
end
