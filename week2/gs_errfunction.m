function E = gs_errfunction(P0, Xobs)
    Xobs = reshape (Xobs,2,[]);
    x = [Xobs(:, 1:size(Xobs,2)/2); ones(1,size(Xobs,2)/2)];
    xp = [Xobs(:, size(Xobs,2)/2+1:end); ones(1,size(Xobs,2)/2)];
    H  = reshape(P0(1:9),3,3);
    
    d1 = computeEuclideanDistance (x,H\xp);
    d2 = computeEuclideanDistance (xp,H*x);
    E = [d1, d2];
end

function dist = computeEuclideanDistance (x1,x2)

    dist = (x1(1,:)./x1(3,:)-x2(1,:)./x2(3,:)).^2+(x1(2,:)./x1(3,:)-x2(2,:)./x2(3,:)).^2;

end