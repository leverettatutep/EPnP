function v = reorderV(V1,NumC)
     if NumC == 3
         Cc = reshape(V1,[3,3]);
         v = reshape(Cc',[9,1]);
     else
         Cc = reshape(V1,[4,3]);
         v = reshape(Cc',[12,1]);
     end
end