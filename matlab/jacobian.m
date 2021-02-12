function [outputArg1,outputArg2] = jacobian(u2,A2,alphas)
    outputArg1 = inputArg1;
    outputArg2 = inputArg2;
end


function dc = dcdxyzj(xyz,j)
    dc = zeros(3,4);
    dc(xyz,j) = 1;
end