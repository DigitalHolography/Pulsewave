function avgA = maskedAverage(A,d,mask)
    % Returns the average of A considering only a masked region
    B = ones(d); % Averaging conv2 kernel
    SumA = conv2(A.*mask,B,"same"); % sum of masked A pixels in range in each pix
    nbnz = conv2(mask,B,"same"); % nb of non null pixels in range

    avgA = SumA./nbnz;

end