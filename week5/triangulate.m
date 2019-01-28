%TRIANGULATE Triangulates a match x1 <--> x2 given the camera matrices.

function X = triangulate(x1, x2, P1, P2, imsize)

% Precondition.
if nargin > 2
  H = [2/imsize(1) 0 -1
       0 2/imsize(2) -1
       0 0            1];
  P1 = H * P1;
  P2 = H * P2;
  x1 = euclid(H * homog(x1));
  x2 = euclid(H * homog(x2));
end

% Create the design matrix A.
A = [ x1(1) * P1(3,:) - P1(1,:)
      x1(2) * P1(3,:) - P1(2,:)
      x2(1) * P2(3,:) - P2(1,:)
      x2(2) * P2(3,:) - P2(2,:) ];

% Solve AX = 0.
[U,D,V] = svd(A,0);
X = V(:,end);

end

