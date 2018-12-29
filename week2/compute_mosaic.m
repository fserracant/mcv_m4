function I = compute_mosaic(imargb, imbrgb, imcrgb, threshold, norm, showfigures)
    if (norm)    
        ima = sum(double(imargb), 3) / 3 / 255;
        imb = sum(double(imbrgb), 3) / 3 / 255;
        imc = sum(double(imcrgb), 3) / 3 / 255;
    else
        ima = imargb/255;
        imb = imbrgb/255;
        imc = imcrgb/255;
    end
    
    %% Compute SIFT keypoints
    [points_a, desc_a] = sift(ima, 'Threshold', threshold, 'Verbosity', 1);
    [points_b, desc_b] = sift(imb, 'Threshold', threshold, 'Verbosity', 1);
    [points_c, desc_c] = sift(imc, 'Threshold', threshold, 'Verbosity', 1);
    
    if (showfigures)
        figure;
        imshow(imargb);%image(imargb)
        hold on;
        plot(points_a(1,:), points_a(2,:),'+y');
        figure;
        imshow(imbrgb);%image(imbrgb);
        hold on;
        plot(points_b(1,:), points_b(2,:),'+y');
        figure;
        imshow(imcrgb);%image(imcrgb);
        hold on;
        plot(points_c(1,:), points_c(2,:),'+y');
    end

    matches_ab = siftmatch(desc_a, desc_b);
    matches_bc = siftmatch(desc_b, desc_c);
    if (showfigures)
        figure;
        plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');
        figure;
        plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');
    end
    

    th = 3;
    xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
    xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
    [Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000);
    
    if (showfigures)
        figure;
        plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
            matches_ab(:,inliers_ab), 'Stacking', 'v');

        vgg_gui_H(imargb, imbrgb, Hab);
    end
    
    xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
    xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
    [Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 
    
    if (showfigures)
        figure;
        plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
            matches_ab(:,inliers_ab), 'Stacking', 'v');

        vgg_gui_H(imargb, imbrgb, Hab);
    end

    corners = [-400 1200 -100 650];
    iwb = apply_H_v2(imbrgb, eye(3), corners);   % ToDo: complete the call to the function
    iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
    iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % ToDo: complete the call to the function
    
    I = max(iwc, max(iwb, iwa));
end
