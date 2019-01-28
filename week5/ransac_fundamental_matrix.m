function [F, idx_inliers] = ransac_fundamental_matrix(x1, x2, th)

[Ncoords, Npoints] = size(x1);

it = 0;
max_it = 1000;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it
    
    points = randomsample(Npoints, 8);
    F = fundamental_matrix(x1(:,points), x2(:,points));
    inliers = compute_inliers(F, x1, x2, th);
      
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(best_inliers) / Npoints;
    pNoOutliers = 1 - fracinliers^8;
    pNoOutliers = max(eps, pNoOutliers);   % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers); % avoid division by 0
    max_it = min(max_it, log(1-p)/log(pNoOutliers));
    
    it = it + 1;

end

% compute F from all the inliers
F = fundamental_matrix(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;


%--------------------------------------------------------------------------
% Function to evaluate the first order approximation of the geometric error
% (Sampson distance) of the fit of a fundamental matrix with respect to a
% set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
% 'Multiple View Geometry in Computer Vision', page 270 (or page 287, 
% eq. (11.9) in the 2nd edition).

function idx_inliers = compute_inliers(F, x1, x2, th)

    x2tFx1 = zeros(1,length(x1));
	for n = 1:length(x1)
	    x2tFx1(n) = x2(:,n)'*F*x1(:,n);
	end
	
	Fx1 = F*x1;
	Ftx2 = F'*x2;     
	
	% evaluate distances
	d =  x2tFx1.^2 ./ ...
	     (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
	
	idx_inliers = find(abs(d) < th);   
   
        
 
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat

    