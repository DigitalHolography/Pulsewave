function mask = magicwand(inputMask, tolerance, conn, numLargestFeatsToKeep)%
% ENTRY :
% im = maskArtery (binarized)
% conn = 8 by default (fonction bwlabel)
% magic wand selects the 'numLargestFeatsToKeep' most connected, largest
% patterns in the inputMask

A = logical(inputMask);
B = double(bwlabel(A,conn));

mask = zeros(size(inputMask));
area = 1:max(B(:));
for ii = 1:max(B(:))
    idxList = find(B==ii);
    area(ii) = length(idxList);
end

figure(111111)
title('histogram magic wand')
plot(area)

for jj = 1:numLargestFeatsToKeep
    [~,idx] = sort(area,'descend');
    mask = mask + (B == idx(jj));
end

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