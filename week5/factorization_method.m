function [Pproj, Xproj] = factorization_method(x, P, lambda_init)
%FACTORIZATION_METHOD computes a
% projective reconstruction with the factorization method of Sturm and
% Triggs '1996
% This function returns an estimate of:
%       Pproj: 3*Ncam x 4 matrix containing the camera matrices
%       Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
% 
% As a convergence criterion you may compute the Euclidean
% distance (d) between data points and projected points in both images 
% and stop when (abs(d - d_old)/d) < 0.1 where d_old is the distance
% in the previous iteration.

[x1,T1]=normalise2dpts(x1);
[x2,T2]=normalise2dpts(x2);

x = [x1;x2];
if strcmpi(lambda_init, 'ones')
    lambda = ones(2,length(x));
elseif strcmpi(lambda_init, 'StrumAndTriggs')
    
end

%x = normalize_data(x);

d = 1000;
converged = false;
while not(converged)
    
    lambda = rescaleLambda(lambda);
    M = zeros(6, length(x2));
    M(1:3,:) = x(1:3,:).*lambda(1,:);
    M(4:6,:) = x(4:6,:).*lambda(2,:);
    
    [U, D, V] = svd(M);

    D_4 = D;
    D_4(D_4 < D_4(4,4)) = 0;

    P_M = U * D_4;
    X_M = V';

    d_old = d;
    d = sum(sqrt( sum((x - (P_M * X_M)) .^ 2)))
    converged = (abs(d - d_old) / d) < 0.1;
    if not(converged)
        xp = P_M*X_M;
        lambda(1,:) = xp(3,:);
        lambda(2,:) = xp(6,:);
    end
end

Pproj = P_M;
Xproj = X_M;

Pproj(1:3,:) = T1 \ Pproj(1:3,:);
Pproj(4:6,:) = T2 \ Pproj(4:6,:);

end

% function x_norm = normalize_data(x)
% 
% x_norm = [];
% for i = 1:3:size(x,1)
%     x_i = x(i:i+2,:);
%     centroid = mean(x_i,2);
%     distances = sqrt( sum((x_i - centroid) .^ 2 ));
%     s_i = sqrt(2) / mean(distances);
%     H_i = [s_i 0 -s_i*centroid(1); 0 s_i -s_i*centroid(2); 0 0 1];
%     x_i_norm = H_i * x_i;
%     
%     x_norm = [x_norm; x_i_norm];
% end
% 
% end


function newLambda = rescaleLambda(lambda)


converged = false;
for i=1:5
    oldLambda = lambda;
    lambda = lambda./vecnorm(lambda,1,1);
    converged = (abs(sum(lambda(:) - oldLambda(:))) / sum(lambda(:))) < 0.1;
    if converged
            break;
    end
    
    oldLambda = lambda;
    lambda = lambda./vecnorm(lambda,1,2);
    converged = (abs(sum(lambda(:) - oldLambda(:))) / sum(lambda(:))) < 0.1;
    if converged
            break;
    end
    
end

newLambda = lambda;

end
