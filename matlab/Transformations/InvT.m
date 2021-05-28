function [tinv,scale] = invT(t)
%     t(1:3,4) = 0;
    scale = 1/nthroot(det(t),3);
    r = t;
    r(1:3,4) = 0;
    tinv = r';
    tinv(:,4) = - r' * t(:,4);
    tinv = tinv*scale^2;
    tinv(4,4) = 1;
end