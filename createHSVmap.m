function [hue,sat,val,cmap] = createHSVmap(Im,mask,hue_val_min,hue_val_max)

%Im grayscale image with value between 0 and 1

hue = (Im*(hue_val_max-hue_val_min) + hue_val_min).*mask;
sat = (1.0-0.5*Im).*mask;
val = Im;
tolVal = [0.02, 0.98];
val = mat2gray(imadjust(val, stretchlim(val, tolVal)));

ARI_x = linspace(0,1,256);
ARI_h = (ARI_x*(hue_val_max-hue_val_min) + hue_val_min);
ARI_s = ones(1,256);
ARI_v = ARI_x;
cmap = squeeze(hsv2rgb(ARI_h,ARI_s,ARI_v));

end

%sat = (1.0 - abs(0.5 * mat2gray(v))) .* double(or(maskArtery,maskVein)) ;