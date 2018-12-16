function J = apply_H(I, H, crop_to_fit)
    switch nargin
      case 2
          crop_to_fit = false;
      case 3
      otherwise
        error('Incorrect input parameters')
    end
    % compute the size of the resulting image after H is applied
    % by applying H to the 4 corners of I and geting height and
    % width from transformed corners
    SI = size(I);
    C = [1 1 1; 1 SI(2) 1; SI(1) 1 1; SI(1) SI(2) 1];
    C_H = (H * C')';
    C_H_x = C_H(:,1) ./ C_H(:,3);
    C_H_y = C_H(:,2) ./ C_H(:,3);
    min_x = floor(min(C_H_x));
    max_x = ceil(max(C_H_x));
    min_y = floor(min(C_H_y));
    max_y = ceil(max(C_H_y));
    if crop_to_fit
        cols = max([max_y,0]) - min([min_y,0])+1;
        rows = max([max_x,0]) - min([min_x,0])+1;
    else

        cols = max_y - min_y+1;
        rows = max_x - min_x+1;
    end
    
    % empty destination image
    J = zeros(rows, cols, 3);
    
    % prepare a matrix with all points in projective coords (i.e. each
    % row is [x y 1])
    [Y,X] = meshgrid(min_y:max_y, min_x:max_x);
    % TODO: correct placing of homography and apply t_x and t_y
    nx = max_x - min_x+1;
    ny = max_y - min_y+1;
    X = X(:);
    Y = Y(:);
    P = [ X Y repmat(1,nx*ny,1)];
    
    % apply inverse homography to all points and separate transformed
    % x and y coords for interp2
    E = H\P';
    E_y = E(1,:); 
    E_x = E(2,:);
    
    % apply interp2 to each channel in original image
    ID = im2double(I);
    if crop_to_fit
        J(max(1,min_x):max(1,min_x)+nx-1,max(1,min_y):max(1,min_y)+ny-1,1) = im2uint8( reshape( interp2(ID(:,:,1), E_x, E_y), [nx,ny]) );
        J(max(1,min_x):max(1,min_x)+nx-1,max(1,min_y):max(1,min_y)+ny-1,2) = im2uint8( reshape( interp2(ID(:,:,2), E_x, E_y), [nx,ny]) );
        J(max(1,min_x):max(1,min_x)+nx-1,max(1,min_y):max(1,min_y)+ny-1,3) = im2uint8( reshape( interp2(ID(:,:,3), E_x, E_y), [nx,ny]) );
    else
        J(1:nx,1:ny,1) = im2uint8( reshape( interp2(ID(:,:,1), E_x, E_y), [nx,ny]) );
        J(1:nx,1:ny,2) = im2uint8( reshape( interp2(ID(:,:,2), E_x, E_y), [nx,ny]) );
        J(1:nx,1:ny,3) = im2uint8( reshape( interp2(ID(:,:,3), E_x, E_y), [nx,ny]) );
    end

    imshow(uint8(J));
    
    
end