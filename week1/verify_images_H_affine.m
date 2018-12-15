function areEqual = verify_images_H_affine(I, H_dec, I2)
% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before

% Compute I2 transformed by the decomposed representation of the affine
% homography H_dec
I2_dec = apply_H(I, H_dec);
% Aside of visually verifying it, we can do it numerically.
% WARNING: maybe they are not exactly equal due to rounding errors. If
% so(after visual check), compute difference against a small threshold
% (1e-7 or similar).
areEqual = isequal(I2, I2_dec);