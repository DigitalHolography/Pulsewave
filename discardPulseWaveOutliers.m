function [idxOut,pulseWaveOutput] = discardPulseWaveOutliers(pulseWave,filterSize)

[~,outliers] = rmoutliers(pulseWave,'movmedian',filterSize);
% [~,outliers] = rmoutliers(pulseWave,'movmean',filterSize);
% [~,outliers] = rmoutliers(pulseWave,'grubbs');
idxOut = find(outliers);
if isempty(idxOut)
else
    for tt = 1:length(idxOut)
        if(idxOut(tt) == 1) % first point
            pulseWave(idxOut(tt)) = mean(pulseWave(idxOut(tt):(idxOut(tt)+1)));
        elseif(idxOut(tt) == length(pulseWave))% endpoint
            pulseWave(idxOut(tt)) = mean(pulseWave((idxOut(tt)-1):idxOut(tt)));
        elseif(idxOut(tt) == length(pulseWave)-1)
            pulseWave(idxOut(tt)) = mean(pulseWave((idxOut(tt)-2):idxOut(tt)));
        elseif(idxOut(tt) == 2)
            pulseWave(idxOut(tt)) = mean(pulseWave((idxOut(tt)-1):(idxOut(tt)+1)));
        else
            pulseWave(idxOut(tt)) = mean([ ...
                pulseWave(idxOut(tt)-2) ...
                pulseWave(idxOut(tt)-1) ...
                pulseWave(idxOut(tt)+1) ...
                pulseWave(idxOut(tt)+2)]);
        end
    end
    % averaging may be counterproductive
end
idxOut;
pulseWaveOutput = pulseWave;
end