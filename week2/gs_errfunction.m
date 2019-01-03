function E = gs_errfunction(P0, Xobs)
    Xobs = reshape (Xobs,2,[]);
    x = [Xobs(:, 1:size(Xobs,2)/2); ones(1,size(Xobs,2)/2)];
    xp = [Xobs(:, size(Xobs,2)/2+1:end); ones(1,size(Xobs,2)/2)];
    H  = reshape(P0(1:9),3,3);
    xhat = [reshape(P0(9+1:end),2,[]);ones(1,size(Xobs,2)/2)];
    xhatp = H*xhat;
     
    d1 = computeEuclideanDistance (x,xhat);
    d2 = computeEuclideanDistance (xp,xhatp);
    E = [d1, d2];
end

function dist = computeEuclideanDistance (x1,x2)
    dist = (x1(1,:)./x1(3,:)-x2(1,:)./x2(3,:)).^2+(x1(2,:)./x1(3,:)-x2(2,:)./x2(3,:)).^2;
end