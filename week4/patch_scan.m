function [err, minarg] = patch_scan(patch, image)
  %% Look for a patch in an image (y-axis only). 
  % Returns difference on each test
  err = [];
  w = size(patch);
  [H,W,D] = size(image);
  
%   for n=1:20:(W - w(2))
%     candidate = image(:, n:(n+w(2)-1));
%     err = [err patch_diff(patch, candidate)];
%   end
  
  p1 = rgb2gray(patch);
  p2 = rgb2gray(image);
  err = -abs(conv2(p2, p1, 'valid'));

  minarg = find(err==min(err));
end
