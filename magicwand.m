function finalMask = magicwand(inputMask,refImage, tolerance, conn, numLargestFeatsToKeep)%
% ENTRY :
% im = maskArtery (binarized)
% conn = 8 by default (fonction bwlabel)
% magic wand selects the 'numLargestFeatsToKeep' most connected, largest
% patterns in the inputMask
A = logical(inputMask);

blurred_mask = imgaussfilt(double(refImage),0.03*size(refImage,1),'Padding',0);
% 
figure(1111)
imagesc(blurred_mask)

[centroidY, centroidX] = find(blurred_mask == max(blurred_mask,[],'all'));

radius = 0.13*size(A,1);
mask = false(size(A));

[rows, cols] = size(mask);
[x, y] = meshgrid(1:cols, 1:rows);
mask = sqrt((x - centroidX).^2 + (y - centroidY).^2) <= radius;
outputImage = A | mask;
figure, imagesc(outputImage)

labeledImage = double(bwlabel(outputImage, conn));
finalMask = zeros(size(outputImage));
area = 1:max(labeledImage(:));
for ii = 1:max(labeledImage(:))
    pixelIndices = find(labeledImage == ii);
    area(ii) = length(pixelIndices);
end
for jj = 1:numLargestFeatsToKeep
    [~, idx] = sort(area, 'descend');
    finalMask = finalMask + (labeledImage == idx(jj));
end
finalMask = finalMask - mask;
finalMask = finalMask + (A & mask);




% 
% A = logical(inputMask);
% B = double(bwlabel(A,conn));
% 
% mask = zeros(size(inputMask));
% area = 1:max(B(:));
% for ii = 1:max(B(:))
%     idxList = find(B==ii);
%     area(ii) = length(idxList);
% end
% 
% figure(111111)
% title('histogram magic wand')
% plot(area)
% 
% for jj = 1:numLargestFeatsToKeep
%     [~,idx] = sort(area,'descend');
%     mask = mask + (B == idx(jj));
% end

%% Method with flip
% maxiB = max(B(:));
% B = B ./ maxiB;
% B = B > tolerance;
% mask1 = xor(A,B);
% 
% AA = flip(A,2);
% BB = double(bwlabel(AA,conn));
% % figure(777)
% % imagesc(BB)
% % colorbar
% maxiBB = max(BB(:));
% BB = BB/maxiBB;
% BB = BB>tolerance;
% mask2 = xor(AA,BB);
% 
% mask = mask1.*flip(mask2,2);
% % figure(888)
% % imagesc(mask)
% % colorbar
% end

end