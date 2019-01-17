function X = triangulate(x1, x2, P1, P2, imsize)
    %TRIANGULATE performs a triangulation with the homogeneous algebraic method (DLT)
    %   The entries are (x1, x2, P1, P2, imsize), where:
    %           - x1, and x2 are the Euclidean coordinates of two matching 
    %             points in two different images.
    %           - P1 and P2 are the two camera matrices
    %           - imsize is a two-dimensional vector with the image size


    H = [2/imsize(1) 0 -1; 0 2/imsize(2) -1; 0 0 1];
    
    x1 = H * [x1; 1];
    x2 = H * [x2; 1];
    P1 = H * P1;
    P2 = H * P2;
    
    A = [ x1(1) * P1(3,:) - P1(1,:) ; ...
          x1(2) * P1(3,:) - P1(2,:) ; ...
          x2(1) * P2(3,:) - P2(1,:) ; ...
          x2(2) * P2(3,:) - P2(2,:) ]
      
    [~, ~, V] = svd(A);
    X = V(:, end)
    
    X = X ./ X(4);
end

