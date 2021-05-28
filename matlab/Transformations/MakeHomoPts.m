function colPts = MakeHomoPts(colVecs)
    colPts = [colVecs; ones(1,size(colVecs,2))];
end