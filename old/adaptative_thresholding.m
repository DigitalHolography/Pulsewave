function [mask_bin] = adaptative_thresholding(img, filtered_img, alpha, nb_div)
    %UNTITLED2 Summary of this function goes here
    
    img_eq = img - min(img,[],"all");
    img_eq = img_eq / max(img_eq,[],"all");
    seg_choro = img_eq > 0.95;
    choro_region = bwareafilt(seg_choro,1);
    stats = regionprops(choro_region,"Centroid","Area");

    radius = stats.Area^0.5/pi;
    cercle_mask = false(size(filtered_img));
    [rows, cols] = size(cercle_mask);
    [x, y] = meshgrid(1:cols, 1:rows);
    circle = sqrt((x - stats.Centroid(1)).^2 + (y - stats.Centroid(2)).^2) <= radius;

    high_seg = filtered_img > (0.5*max(filtered_img,[],'all'));
    high_seg_one_region = bwareafilt(high_seg,10);

    init_mask = init_region_growing(filtered_img, nb_div);

    init_mask = init_mask | high_seg_one_region;


    mask_bin = region_growing_for_vessel(filtered_img, init_mask, circle, alpha);
end