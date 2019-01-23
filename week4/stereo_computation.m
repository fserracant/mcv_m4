
function D = stereo_computation(left, right, min_d, max_d, winsize, matching_cost)
%STEREO_COMPUTATION% computes the disparity between a pair of
% rectified images using a local method based on a matching cost
%   between two local windows.
%
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)

if ndims(left) == 3
  left = im2double(rgb2gray(left)) * 255.0;
end
if ndims(right) == 3
  right = im2double(rgb2gray(right)) * 255.0;
end

D = zeros(size(left));

if mod(winsize, 2) == 0
  winsize = winsize - 1;
end
padding = (winsize - 1) / 2;

[left_rows, left_cols] = size(left);

left = padarray(left, [padding, padding], 'symmetric');%0);
right = padarray(right, [padding, padding], 'symmetric');%0);
% Mask to avoid using padded values (set to 0) in matching costs
%     mask_p = ones(size(left));
%     mask_p = padarray(mask_p, [padding, padding], 0);

gamma_c = 5.0;
gamma_p = (winsize -1) / 2;%17.5;  % Window's radius
tau_col = 40.0;

if strcmpi(matching_cost, 'BW')

% Based off the IPOL paper
win_r = (winsize-1) / 2;
[X,Y] = meshgrid(-win_r:win_r, -win_r:win_r);
w_pos = exp(-2.0 * sqrt(X(:).^2 + Y(:).^2) / gamma_p);



end

for row = 1+padding:left_rows+padding
  for col = 1+padding:left_cols+padding
    % Define starting and finishing disparity (search range)
    start_disp = max(1+padding, col - max_d);
    end_disp = min(left_cols + padding, col + max_d);
    % Mask evaluated at current window
    left_window = left(row-padding:row+padding, col-padding:col+padding);
    weights = zeros(size(left_window));
    weights(:) = 1/(numel(left_window));
    
    if strcmpi(matching_cost, 'SSD')  % Sum of squared differences
      bestSSD = Inf;
      for d = start_disp:end_disp
        right_window = right(row-padding:row+padding,d-padding:d+padding);
        SSD = sum(sum(weights .* abs(left_window - right_window).^2));
        if SSD < bestSSD
          bestSSD = SSD;
          best_idx = d;
        end
      end
      
    elseif strcmpi(matching_cost, 'SAD')
      bestSAD = Inf;
      for d = start_disp:end_disp
        right_window = right(row-padding:row+padding,d-padding:d+padding);
        SAD = sum(sum(weights .* abs(left_window - right_window)));
        if SAD < bestSAD
          bestSAD = SAD;
          best_idx = d;
        end
      end
    elseif strcmpi(matching_cost, 'NCC')
      bestNCC = -Inf;
      for d = start_disp:end_disp
        right_window = right(row-padding:row+padding,d-padding:d+padding);
        sum_left = sum(weights(:) .* left_window(:));
        sum_right = sum(weights(:) .* right_window(:));
        sigma_left = sqrt(sum(weights(:) .* (left_window(:) - sum_left).^2));
        sigma_right = sqrt(sum(weights(:) .* (right_window(:) - sum_right).^2));
        NCC= sum( weights(:).*(left_window(:)-sum_left).*(right_window(:)-sum_right) )...
          /(sigma_left*sigma_right);
        
        if NCC > bestNCC
          bestNCC = NCC;
          best_idx = d;
        end
      end
    end
      D(row-padding, col-padding) = abs(col - best_idx);
      
    end
  end
 
  D(D < min_d) = min_d;
  D = D - min(D(:));
  D(D > max_d) = max_d;
  D = uint8(255 * mat2gray(D));
  
end

