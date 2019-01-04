function [logo_proj, H2] = insertLogo(newLogo, oldLogo, H, corners)
% Inserts 'newLogo' where 'oldLogo' was in a destination image with
% 'dstcorners'.
%
% Inputs
%   - 'newLogo':        logo to replace 'oldLogo'(i.e.: the one to be inserted).
%   - 'oldLogo':        logo that will be replaced.
%   - 'H':              homography that related the oldLogo & the
%                       destination image.
%   - 'corners':        corners of the destination image where the logo
%                       will be inserted.
%
% Outputs
%   - 'logo_proj':      logo projected via H + scaling transformation.
%   - 'H2':             transformation that maps 'newLogo' to 'logo_proj'.

[dst_height, dst_width, ~] = size(newLogo);
[src_height, src_width, ~] = size(oldLogo);

% Define scale transformation
S = [src_width/dst_width 0 0; 0 src_height/dst_height 0; 0 0 1];

% Compute total transformation as H * S
H2 = H * S;

% Apply transformation to target logo
logo_proj = apply_H_v2(newLogo, H2, corners);
