function tinv = invT(t)
    tinv = t';
    tinv(4,1:3) = 0;
    temp = tinv * t;
    tinv(:,4) = -tinv * t(:,4);
    tinv(4,4) = 1;
end