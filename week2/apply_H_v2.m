function I2=apply_H_v2(I, H, corners)

% Apply the transformation H to the image I
%
% INPUTS:
%   I: image to be transformed (color or gray scale image)
%   H: 3x3 matrix that specifies the desired transformation
%   corners: vector of output image corners
%
% OUTOUT:
%   I2: transformed image
%

% input image size
[nrows, ncols, nchan] = size(I);
 
xmin = corners(1);
xmax = corners(2);
ymin = corners(3);
ymax = corners(4);

% create matrices of homogeneous coordinates
[X,Y] = meshgrid(xmin:xmax,ymin:ymax);
Hncols = xmax - xmin + 1;
Hnrows = ymax - ymin + 1;
Z = ones(Hnrows,Hncols);


% create a 3x(Hnrows*Hncols) matrix in order to transform all the 
% coordinate points with matrix operations (thus avoiding the use of a for loop)
Xs = reshape(X,1,Hnrows*Hncols);
Ys = reshape(Y,1,Hnrows*Hncols);
Zs = reshape(Z,1,Hnrows*Hncols);
XYZs = [Xs; Ys; Zs];

% transform image
Hi = inv(H);
HiXYZs = Hi*XYZs;
HX = reshape(HiXYZs(1,:),Hnrows,Hncols);
HY = reshape(HiXYZs(2,:),Hnrows,Hncols);
HZ = reshape(HiXYZs(3,:),Hnrows,Hncols);
HX = HX./HZ;
HY = HY./HZ;

I2 = zeros(Hnrows, Hncols, nchan);
for c=1:nchan,
    I2(:,:,c) = interp2(double(I(:,:,c)),HX,HY,'linear');
end

I2 = uint8(I2);

