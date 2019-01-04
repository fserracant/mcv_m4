function [H, dstcorners] = detectLogos(imdst, imlogo, sift_thr)  
% Detects logo 'im_logo' in image 'im_src' by computing sift matches
% between both images and then estimating an homography via the robust DLT
% algorithm.
%
% Inputs
%   - 'imdst':              destination image where the logo appears somewhere.
%   - 'imlogo':             target logo to detect in the image 'imsrc'.
%   - 'sift_thr':           threshold to filter SIFT matches.
%
% Outputs
%   - 'H':                  estimated homography between 'imsrc' & 'imlogo'.
%   - 'dstcorners':         corners in the logo transformed by 'H'.
%
% About selecting the SIFT threshold:
%
% We need to play with different threshold values in the building case.
% The problem is that the UPF logo on the building facade is not exactly 
% the same logo.
% 0.1  finds no matches at all
% 0.01 finds good points but a bad homography is produced because of logos
%         mismatch
% 0.05 finds some wrong matches but the resulting homography is more
%         accurate to find the corners of the logo

v = 0;  % Set verbosity for SIFT.
%% Compute SIFT descriptors
[points_src, desc_src] = sift(imlogo, 'Threshold', sift_thr, 'Verbosity', v);
[points_dst, desc_dst] = sift(imdst, 'Threshold', sift_thr, 'Verbosity', v);

%% Compute SIFT matches
matches = siftmatch(desc_src, desc_dst);

%% Compute homography via DLT algorithm
ransac_th = 3;
xab_a = [points_src(1:2, matches(1,:)); ones(1, length(matches))];
xab_b = [points_dst(1:2, matches(2,:)); ones(1, length(matches))];
[H, ~] = ransac_homography_adaptive_loop(xab_a, xab_b, ransac_th, 1000);

%% Compute new corners for the transformed logo (via 'H')
[rows, cols, ~] = size(imlogo);
% outter border of the logo
border = 5;
srccorners = [border border 1; border rows-border 1; cols-border border 1; cols-border rows-border 1];
dstcorners = H * srccorners';
dstcorners = dstcorners ./ dstcorners(3,:);
