function avgA = maskedAverage(A, d, maskzone, maskapply)

% Returns the average of A considering only a masked region maskzone
% and apply it on the maskapply

B = ones(d); % Averaging conv2 kernel
SumA = conv2(A .* maskzone, B, "same"); % sum of masked A pixels in range in each pix
nbnz = conv2(maskzone, B, "same"); % nb of non null pixels in range
nbnz((nbnz) == 0) = 1;
avgA = maskapply .* SumA ./ nbnz + ~maskapply .* A;

end
