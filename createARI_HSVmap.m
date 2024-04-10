function [hue,sat,val] = createARI_HSVmap(ARImap,Im,mask,ToolBox)


ARImap = imgaussfilt(ARImap,1.3);
hue = ((ToolBox.ARI_hue_max-ToolBox.ARI_hue_min)./(1+exp(-(ARImap-ToolBox.ARI_inflexion_point_hue).*ToolBox.ARI_slope_hue)) + ToolBox.ARI_hue_min).*mask;


sat = ones(size(Im,1),size(Im,2)).*mask ; %(1.0-0.5*Im).*mask;


val_artery = (((ToolBox.ARI_val_max-ToolBox.ARI_val_min)./(1+exp(-(ARImap-ToolBox.ARI_inflexion_point_val).*ToolBox.ARI_slope_val)) + ToolBox.ARI_val_min)).*mask;


val_bckg = Im;
tolVal = [0.02, 0.98];
val_bckg = mat2gray(imadjust(val_bckg, stretchlim(val_bckg, tolVal))).*(~mask);
val = val_bckg + val_artery;



end