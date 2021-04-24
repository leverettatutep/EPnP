function Cc_ = reshapeVector(EigenVector)
    numRows = size(EigenVector,1);
    if numRows == 12
        Cc_ = reshape(EigenVector,[3,4])'; %Cc row 1 is CP 1, Cc row 2 is CP 2 
    else
        Cc_ = reshape(EigenVector,[3,3])';
    end
end

