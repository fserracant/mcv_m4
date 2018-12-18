function areEqual = verify_product_H_affine(H, H_composed)
% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

tolerance = 1e-6;
areEqual = sum(H(:) - H_composed(:) < tolerance)==numel(H);
