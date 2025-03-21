function [diameter] = findPapilla(M0img, xybary, maskVessels)
%Returns the diameter of the papilla measured

[numX, numY] = size(M0img);

M0hide = maskedAverage(M0img, 30, ~maskVessels);

M0hide = M0hide(xybary(1) + (-round(numX / 4):round(numX / 4)), xybary(2) + (-round(numY / 4):round(numY / 4)));
[Gx, Gy] = gradient(M0hide);
Sob = imadjust(abs(Gx) + abs(Gy));
figure(34); imshow(Sob);
% binaryImg = imbinarize(Sob, 'global');

[centers, radii] = imfindcircles(imadjust(M0hide), [80, 200], ...
    'ObjectPolarity', 'dark', 'Sensitivity', 0.95);
figure(45), imshow(imadjust(M0hide)); title('M0hide'); % Show original image
viscircles(centers, radii, 'EdgeColor', 'b'); % Overlay detected circles
diameter = radii(1) * 2;
end
