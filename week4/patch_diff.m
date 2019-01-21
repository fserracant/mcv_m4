function [error] = patch_diff(patch1, patch2)
  % Euclidean distance between patches
  diff = sqrt(sum(sum((patch1-patch2).^2)));
%   crr = xcorr2(rgb2gray(patch1), rgb2gray(patch1))
%   [ssr,snd] = max(crr(:))
%   [ij,ji] = ind2sub(size(crr),snd)
%   diff = sum(diff);
  % 
%   p1 = rgb2gray(patch1);
%   p2 = rgb2gray(patch2);
%   error = -sum(conv2(patch2, patch1, 'same'));
end