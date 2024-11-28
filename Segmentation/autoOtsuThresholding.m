function [maskArtery, maskVein, maskChoroid] = autoOtsuThresholding(image, mask, classes, name, ToolBox)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% classes = [0 0 1 1] a 0 and 1 matrix where 1 is the class selected and 0 the rejected

cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;
cChoroid = [0 179 0] / 255;

numClasses = size(classes, 1);
color = zeros(3, numClasses);
for i = 1:numClasses
    switch classes(i)
        case 1
            color(:, i) = cArtery;
        case -1
            color(:, i) = cVein;
        case 0
            color(:, i) = [1 1 1];
        case 2
            color(:, i) = cChoroid;
    end
end

level = multithresh(image(mask), numClasses - 1);
graphThreshHistogram(image, level, mask, color', name, ToolBox)
level = [-1 level];
quantizedImage = imquantize(image - 2 * ~mask, level);
imwrite(rescale(quantizedImage), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s_Quantize.png", ToolBox.main_foldername, name)))

maskArtery = zeros(size(image));
maskVein = zeros(size(image));
maskChoroid = zeros(size(image));

for i = 1:numClasses

    if classes(i) == 1
        maskArtery = maskArtery + (quantizedImage == i + 1);
    elseif classes(i) == -1
        maskVein = maskVein + (quantizedImage == i + 1);
    elseif classes(i) == 2
        maskChoroid = maskChoroid + (quantizedImage == i + 1);
    end

end

end
