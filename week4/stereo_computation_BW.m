
function D = stereo_computation_BW(left, right, min_d, max_d, winsize, matching_cost)
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
        left = im2double(rgb2gray(left)) ;
    end
    if ndims(right) == 3
        right = im2double(rgb2gray(right));
    end
    
    D = zeros(size(left));
    
    if mod(winsize, 2) == 0
        winsize = winsize - 1;
    end
    padding = (winsize - 1) / 2;
    
    left = padarray(left, [padding, padding], 0);
    right = padarray(right, [padding, padding], 0);
    [left_rows, left_cols] = size(left);

    gamma_c = 5;
    gamma_p = winsize/2;
    T = .01;

    if strcmpi(matching_cost, 'BW')
        % bilateral weights
        win_center_idx = floor(winsize/2);
        [X,Y] = meshgrid(-win_center_idx:win_center_idx, -win_center_idx:win_center_idx);
        bw_g = sqrt( (X(:)-X(win_center_idx,win_center_idx)).^2 + ...
                 (Y(:)-Y(win_center_idx,win_center_idx)).^2 );
        bw_g = exp(bw_g / (-gamma_p)); 
        
        %figure;
        %surf(X,Y,reshape(bw_g, winsize, winsize));
    end
    
    for row = 1+padding:left_rows-padding
        right_strip = right(row-padding:row+padding,:);
        right_strip_cols = im2col(right_strip,[winsize winsize]);
            
        for col = 1+padding:left_cols-padding
            real_col = col - padding;
            left_window = left(row-padding:row+padding, col-padding:col+padding);

            
            if strcmpi(matching_cost, 'SSD')  % Sum of squared differences
                C = repmat(left_window(:), [1, size(right_strip_cols,2)]) - right_strip_cols;
                C = sum(C.^2);
                best_idx = find(C == min(C), 1, 'First'); 
                
            elseif strcmpi(matching_cost, 'BW')  % Bilateral weights
                % w(p,q) term
                wl = exp( abs(left_window(:) - left_window(win_center_idx,win_center_idx)) ...
                               / ...
                           (gamma_c * -1));
                wl = wl .* bw_g;
                
                % w(pd,qd) term
                center_idx = ceil(size(right_strip_cols,1)/2);
                wr = exp( abs(right_strip_cols - right_strip_cols(center_idx,:)) ...
                               / ...
                           (gamma_c * -1));
                wr = wr .* bw_g;
                
                % c(q,qd) term
                c = abs(right_strip_cols - left_window(:));
                c(c > T) = T;
                
                % Cost function
                C = sum((wr .* wl) .* c) ./ sum(wr .* wl);  
                
                % apply max_d max disparity around current column
                start_disp = max(1, real_col - max_d);
                end_disp = min(size(C,2), real_col + max_d);
                C(1:start_disp) = Inf;
                C(end_disp:end) = Inf;
                
                % index of the column of the min cost
                best_idx = find(C == min(C), 1, 'Last');
            end
            
            D(row-padding, col-padding) = abs(real_col - best_idx);
            
        end
    end

    D = uint8(255 * mat2gray(D));    
end

