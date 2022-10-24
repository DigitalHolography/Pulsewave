function [one_cycle_video,selectedPulseIdx] = create_one_cycle(video, mask, sys_index_list, Ninterp)
%   one_cycle() : identifies pulse cycles and average them to one video
%   sys_index_list : list of systole indexes in video
arguments
    video
    mask
    sys_index_list {mustBeNumeric, mustBeNonnegative} = []
    Ninterp {mustBeNumeric, mustBePositive} = 64
end

disp('interpolate, average, and shift');
disp('arterial pulse waveforms');
disp('over one cardiac cycle...');

M = length(sys_index_list)-1; % M : pulse #
one_cycle_video = zeros(size(video,1), size(video,2), Ninterp, M);
for ii = 1:M % for each detected pulse, loop
    interp_range = linspace(sys_index_list(ii),sys_index_list(ii+1)-1,Ninterp);
    for id_x = 1 : size(video,1)
        for id_y = 1 : size(video,2)
            one_cycle_video(id_x,id_y,:,ii) = interp1((sys_index_list(ii):sys_index_list(ii+1)-1),squeeze(video(id_x, id_y,sys_index_list(ii):sys_index_list(ii+1)-1)),interp_range);
        end
    end
end

%check pulse integrity
    pulseMatrix = zeros(M,Ninterp);
    cleanSignal = squeeze(mean(one_cycle_video(:,:,:,:),[1 2 4]));
    cleanSignal = cleanSignal - mean(cleanSignal(:));
for ii = 1:M % for each detected pulse, loop
    signal = squeeze(mean(one_cycle_video(:,:,:,ii),[1 2]));
    signal = signal - mean(signal(:));
    pulseMatrix(ii,:) = signal;
    %     cleanSignal = smoothdata(signal,'lowess');
    noise = sqrt(abs(abs(signal).^2 - abs(cleanSignal).^2));
    idxOutNoise = find(noise>2*std(noise));
    dataReliabilityIndex(ii) = ceil(100*(1-(length(idxOutNoise)/length(noise))));
    disp(['data reliability index for pulse ', num2str(ii), ' : ', num2str(dataReliabilityIndex(ii)), ' %']);
end

% [U,S,V] = svd(pulseMatrix);
% s1 = S(1,1);
% S = zeros(size(S)); 
% S(1,1) = s1;
% pulseMatrixFiltered = U * S * V';
% [row col] = find(S>.3*max(S(:)))

tmp = zeros(size(video,1), size(video,2), Ninterp);
ctr = 0;
for ii = 1:M
    if (dataReliabilityIndex(ii) < 50)
        tmp = tmp + squeeze(one_cycle_video(:,:,:,ii)); % average selected cycles
        ctr = ctr + 1;
    end
end
% selectedPulseIdx = find(dataReliabilityIndex(ii) < 50);

if(ctr == 0)
    tmp = squeeze(mean(one_cycle_video,4));% average all cycles
else
    tmp = tmp /ctr;
end

% FIXME : create zerolag correlation matrix between all zero-mean pulses,
% then svd, to identify the most correlated pulses

% FIXME : export selected pulse indexes

selectedPulseIdx = 0;



one_cycle_video = tmp; % average all detected cycles
oneP = squeeze(sum(one_cycle_video .* mask,[1 2 ]) / nnz(mask));
[min_val,shift] = min(oneP); % find bottom systole
one_cycle_video = circshift(one_cycle_video, -shift, 3); % shift pulse to start & end with bottom diastole

disp('done.');

end