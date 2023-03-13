function [ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma] = construct_resistivity_index(video, maskArtery)
% https://en.wikipedia.org/wiki/Arterial_resistivity_index
% arterial resistivity

video = double(video);

arterialPulse = squeeze(sum(video .* maskArtery,[1 2]));
arterialPulse = arterialPulse/nnz(maskArtery);

%% avg. arterial resistivity index

for ii=1:size(video,1)
    for kk=1:size(video,2)
        if maskArtery(ii,kk)==1
            video(ii,kk,:) = imgaussfilt(video(ii,kk,:),2);
        end
    end
end

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
% sat = abs(ARImap .* maskArtery);
% sat = satAmp * sat;%imadjust(sat, stretchlim(sat, tolSat));

% hue = 1 * ones(size(video,1), size(video,2)); % pure red color

% val = flat_field_correction(mat2gray(squeeze(mean(video,3))), 0.07*size(video,1), 0);
% val = mat2gray(val);
% tolVal = [0.02, 0.98]; 
% lowhigh = stretchlim(val, tolVal); % adjust contrast a bit 
% val = imadjust(val, stretchlim(val, tolVal));

% ARImapRGB = hsv2rgb(hue, sat, val);


val = flat_field_correction(mat2gray(squeeze(mean(video,3))), 0.07*size(video,1), 0);
val = mat2gray(val);
tolVal = [0.02, 0.98]; 
lowhigh = stretchlim(val, tolVal); % adjust contrast a bit 
val = mat2gray(imadjust(val, stretchlim(val, tolVal)));

Artery_val = abs(ARImap .* maskArtery);
maxArtery = max(Artery_val(:));
minArtery = min(nonzeros(Artery_val));
for ii=1:size(ARImap,1)
    for kk=1:size(ARImap,2)
        if Artery_val(ii,kk) ~= 0
            Artery_val(ii,kk) = (Artery_val(ii,kk)-minArtery)./(maxArtery-minArtery);
        end
    end
end

Artery_val = Artery_val.*maskArtery;

gamma = 0.5;
ARImapRGB = ones(size(ARImap,1), size(ARImap,2), 3);
ARImapRGB(:,:,1) = val - maskArtery.*val + ones(size(ARImap,1), size(ARImap,2)).*maskArtery;
ARImapRGB(:,:,2) = val - maskArtery.*val + Artery_val.^gamma;
ARImapRGB(:,:,3) = val - maskArtery.*val + Artery_val.^gamma;

%% arterial resistivity video RGB
video = mat2gray(video); % not quantitative anymore
sat = satAmp * mat2gray(abs(ARImap)) .* maskArtery;
% sat = satAmp * abs(ARImap .* maskArtery);
hue = 1 * ones(size(video,1), size(video,2)); % pure red color

%% FIXME : TO TEST
% hue = mat2gray(squeeze(mean(v_RMS,3)))*0.18 .* maskArtery;
% sat = 0.75 * double(maskArtery);
% val = squeeze(mean(v_RMS,3));
% val = mat2gray(val);
% tolVal = [0.02, 0.98]; 
% lowhigh = stretchlim(val, tolVal); % adjust contrast a bit 
% val = mat2gray(imadjust(val, stretchlim(val, tolVal)));

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