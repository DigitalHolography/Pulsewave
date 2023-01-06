function [ARImap, ARI, ARImapRGB, ARIvideoRGB] = construct_resistivity_index(video, maskArtery)
% https://en.wikipedia.org/wiki/Arterial_resistivity_index
% arterial resistivity

video = double(video);

arterialPulse = squeeze(sum(video .* maskArtery,[1 2]));
arterialPulse = arterialPulse/nnz(maskArtery);

%% avg. arterial resistivity index

[maxAP,maxAPidx] = max(arterialPulse(:));
[minAP,minAPidx] = min(arterialPulse(:));

ARI = (maxAP - minAP)/maxAP;

%% arterial resistivity map values (no mask)
% FIXME add time blur to video 
ARImap = squeeze( (video(:,:,maxAPidx) - video(:,:,minAPidx)) ./ video(:,:,maxAPidx) );



%%
% composite image parameters
satAmp = 0.75; %MUST be in [0 1] 

%% arterial resistivity map RGB
sat = abs(ARImap .* maskArtery);
sat = satAmp * sat;%imadjust(sat, stretchlim(sat, tolSat));
hue = 1 * ones(size(video,1), size(video,2)); % pure red color
val = flat_field_correction(mat2gray(squeeze(mean(video,3))), 0.07*size(video,1), 0);
val = mat2gray(val);
tolVal = [0.02, 0.98]; 
lowhigh = stretchlim(val, tolVal); % adjust contrast a bit 
val = imadjust(val, stretchlim(val, tolVal));
ARImapRGB = hsv2rgb(hue, sat, val);

%% arterial resistivity video RGB
video = mat2gray(video); % not quantitative anymore
sat = satAmp * abs(ARImap .* maskArtery);
hue = 1 * ones(size(video,1), size(video,2)); % pure red color
% for pp = 1:size(video,3)
%     video(:,:,pp) = flat_field_correction(squeeze(video(:,:,pp)), 0.07*size(video,1), 0);
% end
adjustedVideo = mat2gray(video);
avgAdjustedVideo = squeeze(mean(adjustedVideo,3));
tolVal = [0.1, 0.99]; 
lowhighVal = stretchlim(avgAdjustedVideo, tolVal); % adjust video contrast a bit 
% tolSat = [0.0, 1]; 
% lowhighSat = stretchlim(avgAdjustedVideo .* sat, tolSat); 
for pp = 1:size(video,3)
    adjustedVideo(:,:,pp) = imadjust(squeeze(adjustedVideo(:,:,pp)), lowhighVal);
    adjustedVideo(:,:,pp) = squeeze(adjustedVideo(:,:,pp));
    adjustedSat(:,:,pp) = squeeze(adjustedVideo(:,:,pp) .* sat); 
end
ARIvideoRGB = uint8(zeros(size(video,1), size(video,2), 3, size(video,3)));
% ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
for pp = 1:size(video,3)
    ARIvideoRGB(:,:,:,pp) = im2uint8(hsv2rgb(hue, adjustedSat(:,:,pp), squeeze(adjustedVideo(:,:,pp))));
end

end