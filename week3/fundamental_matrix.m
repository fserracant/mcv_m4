% 

function F = fundamental_matrix(x1, x2)

    % normalisation of points
    [x1_norm, T1] = normalise2dpts(x1)
    [x2_norm, T2] = normalise2dpts(x2);
    
    % construction of W in 8-point algorithm
    u1 = x1_norm(1,:)';
    v1 = x1_norm(2,:)';
    u2 = x2_norm(1,:)';
    v2 = x2_norm(2,:)';
    
    W = [ u1.*u2 v1.*u2 u2 u1.*v2 v1.*v2 v2 u1 v1 ones(length(u1),1) ];
    
    [~, ~, V] = svd(W,0);
    % Reshape last columnt of V into a 3x3 matrix F
    F = reshape(V(:, end), 3, 3)';
    
    % F may be of rank 3, let's force it to be of rank 2
    [U, D, V] = svd(F,0);  
    D(3,3) = 0;
    
    % reconstruct F with rank 2
    F = U * D * V';
    
    % de-normalise F
    F = T2' * F * T1;
    

    