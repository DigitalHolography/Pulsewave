function [] = ArterialResistivityIndex(onePulseVideo,dataCubeM2M0, maskArtery,  ToolBox)

disp('arterial resistivity...');

meanIm = mat2gray(squeeze(mean(dataCubeM2M0,3)));
tolVal = [0.02, 0.98]; 
meanIm = mat2gray(imadjust(meanIm, stretchlim(meanIm, tolVal)));

[ARI, ARImap] = construct_resistivity_index(onePulseVideo, maskArtery);


[hue_ARI,sat_ARI,val_ARI] = createARI_HSVmap(ARImap,meanIm,maskArtery,ToolBox);
% arterial resistivity map RGB
ARImapRGB = hsv2rgb(hue_ARI, sat_ARI, val_ARI);

disp('done.');

% disp('arterial resistivity...');
% [ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma, img_avg] = construct_resistivity_index(onePulseVideo, maskArtery,path);
% ARImap = ARImap.*maskArtery;
% 
% w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ARIvideoRGB.avi')));
% open(w)
% ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
% for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
%     writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
% end
% close(w);
% 
% disp('done.');


%% Display Figure 
figure(70)
imagesc(ARImapRGB) ;
title(strcat('Arterial resistivity. avg. index value : ', num2str(ARI)));
axis image
axis off
set(gca,'LineWidth', 2);
fontsize(gca,12,"points") ;
c = colorbar('southoutside','Ticks',linspace(0,1,11));
c.Label.String = 'Arterial resistivity index';
c.Label.FontSize = 12;

ARI_x = linspace(0,1,256);
ARI_h = (ToolBox.ARI_hue_max-ToolBox.ARI_hue_min)*sigmoid(ARI_x,ToolBox.ARI_inflexion_point_hue,ToolBox.ARI_slope_hue)+ ToolBox.ARI_hue_min;
ARI_s = ones(1,256);
ARI_v = (ToolBox.ARI_val_max-ToolBox.ARI_val_min)*sigmoid(ARI_x,ToolBox.ARI_inflexion_point_val,ToolBox.ARI_slope_val)+ToolBox.ARI_val_min;
cmap = squeeze(hsv2rgb(ARI_h,ARI_s,ARI_v));
colormap(cmap);

    % Save colorbar
colorfig = figure(113);
colorfig.Units = 'normalized';
colormap(cmap)
hCB = colorbar('north');
set(gca,'Visible',false)
set(gca,'LineWidth', 3);
hCB.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca,15,"points") ;

%% Save Figures
print('-f113','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_ARI_map_colorbar.png')));
print('-f70','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_ARI_map.png')));

print('-f70','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_resistivityMap.eps'))) ;

imwrite(ARImapRGB,fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_ARI_map_img.png')),'png');

close all 
