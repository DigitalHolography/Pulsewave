function [maskSection, VesselImageRGB] = createMaskSection(TB, img, r1, r2, xy_barycenter, figname, maskArtery, maskVein, opt)

arguments
    TB
    img
    r1
    r2
    xy_barycenter
    figname
    maskArtery
    maskVein = []
    opt.thin = 0
end

[numX, numY] = size(maskArtery);
x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;

if opt.thin > 0
    maskSection1 = diskMask(numX, numY, r1, r1 + opt.thin, center = [x_c, y_c]);
    maskSection2 = diskMask(numX, numY, r2 - opt.thin, r2, center = [x_c, y_c]);
    maskSection = xor(maskSection1, maskSection2);
else
    maskSection = diskMask(numX, numY, r1, r2, center = [x_c, y_c]);
end

%% Create Colormap Artery/Vein
img = mat2gray(img);
img_RGB = cat(3, img, img, img);
maskSectionRGB = cat(3, maskSection, maskSection, maskSection);

maskArterySection = maskArtery & maskSection;

cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapArterySection = cmapLAB(256, [1 1 1], 0, [1 1 0], 1/3, [1 0 0], 2/3, [0 0 0], 1);

M0_Artery = setcmap(img, maskArtery, cmapArtery);
M0_ArterySection = setcmap(img, maskArterySection, cmapArterySection);

if isempty(maskVein)

    VesselImageRGB = ...
        M0_Artery .* ~maskSection + ...
        M0_ArterySection + ...
        maskSectionRGB .* ~maskArtery + ...
        img_RGB .* ~maskArtery .* ~maskSection;

    if ~isfolder(fullfile(TB.path_png, 'mask'))
        mkdir(fullfile(TB.path_png, 'mask'));
    end

    imwrite(VesselImageRGB, fullfile(TB.path_png, 'mask', sprintf("%s_%s.png", TB.main_foldername, figname)), 'png');

else

    maskVeinSection = maskVein & maskSection;
    maskAV = maskArtery & maskVein;
    maskAVSection = maskAV & maskSection;

    cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
    cmapVeinSection = cmapLAB(256, [1 1 1], 0, [0 1 1], 1/3, [0 0 1], 2/3, [0 0 0], 1);
    cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);
    cmapAVSection = cmapLAB(256, [0 0 0], 0, [0 1 0], 1/3, [1 1 1], 1);

    M0_Vein = setcmap(img, maskVein, cmapVein);
    M0_VeinSection = setcmap(img, maskVeinSection, cmapVeinSection);
    M0_AV = setcmap(img, maskAV, cmapAV);
    M0_AVSection = setcmap(img, maskAVSection, cmapAVSection);

    maskVessel = maskArtery | maskVein;

    VesselImageRGB = ...
        M0_Artery .* ~maskAV .* ~maskSection + ...
        M0_ArterySection .* ~maskAVSection + ...
        M0_Vein .* ~maskAV .* ~maskSection + ...
        M0_VeinSection .* ~maskAVSection + ...
        M0_AV .* ~maskSection + ...
        M0_AVSection + ...
        maskSectionRGB .* ~maskVessel + ...
        img_RGB .* ~maskArtery .* ~maskVein .* ~maskSection;

    if ~isfolder(fullfile(TB.path_png, 'mask'))
        mkdir(fullfile(TB.path_png, 'mask'));
    end
    imwrite(VesselImageRGB, fullfile(TB.path_png, 'mask', sprintf("%s_%s.png", TB.main_foldername, figname)), 'png');
end

end