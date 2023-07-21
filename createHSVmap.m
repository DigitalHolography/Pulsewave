function [hue,sat,val] = createHSVmap(Im,mask,hue_val_min,hue_val_max)

%Im grayscale image with value between 0 and 1

hue = (Im*(hue_val_max-hue_val_min) + hue_val_min).*mask;
sat = (1.0-0.5*Im).*mask;
val = Im;
tolVal = [0.02, 0.98];
val = mat2gray(imadjust(val, stretchlim(val, tolVal)));

end

%sat = (1.0 - abs(0.5 * mat2gray(v))) .* double(or(maskArtery,maskVein)) ;