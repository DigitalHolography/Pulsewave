function [ARI, ARImap] = construct_resistivity_index(video,maskArtery)
% https://en.wikipedia.org/wiki/Arterial_resistivity_index
% arterial resistivity

video = double(video);

arterialPulse = squeeze(sum(video .* maskArtery,[1 2]));
arterialPulse = arterialPulse/nnz(maskArtery);

%% avg. arterial resistivity index
timeBlurWindow = 7;
video_blurred = video;

for ii = 1:size(video, 1)
    for kk = 1:size(video, 2)
        if maskArtery(ii, kk) == 1
            video_blurred(ii, kk, :) = movmean(imgaussfilt(video(ii,kk,:),3), timeBlurWindow, 3);
        end
    end
end
% figure(549)
% plot(arterialPulse,'k','LineWidth',2)
% set(gca,'LineWidth',2)
% axis tight;


[maxAP,maxAPidx] = max(arterialPulse(:));
[minAP,minAPidx] = min(arterialPulse(:));

ARI = (maxAP - minAP)/maxAP;

%% arterial resistivity map values

ARImap = squeeze( (video(:,:,maxAPidx) - video(:,:,minAPidx)) ./ video(:,:,maxAPidx) );
ARImap(ARImap>1)=1;
ARImap = ARImap .* (ARImap.*maskArtery > 0);

end