% 

function F = fundamental_matrix(x1, x2)
    
    x1_norm = x1 ./ x1(3,:);
    x2_norm = x2 ./ x2(3,:);
    
    u1 = x1_norm(1,:)';
    v1 = x1_norm(1,:)';
    u2 = x2_norm(2,:)';
    v2 = x2_norm(2,:)';
    
    W = [ u1.*u2 v1.*u2 u2 u1.*v2 v1.*v2 v2 u1 v1 ones(length(u1),1) ];
    
    [~, ~, Vt] = svd(W);

    V = Vt';
    f = V(:, end);
    % Reshape h into a 3x3 matrix
    F = reshape(f, 3, 3);
    
    [U, D, V] = svd(F);
    
    Drank2 = D;
    Drank2(3,3) = 0.0;
    
    F = U * Drank2 * V'
    

    