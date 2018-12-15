function areEqual = verify_product_H_affine(H, H_dec)
% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

areEqual = isequal(H, H_dec);
% WARNING: maybe they are not exactly equal due to rounding errors. If
% so(after visual check), compute difference against a small threshold
% (1e-7 or similar).