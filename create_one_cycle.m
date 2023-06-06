function [one_cycle_video,selectedPulseIdx, signal_one_cycle] = create_one_cycle(video, mask, sys_index_list, Ninterp)
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

Reliability_index_threshold = 50;

Nx = size(video,1);
Ny = size(video,2);

if M > 0 % we have detected at least two systoles
    parfor ii = 1:M % for each detected pulse, loop
        interp_range = linspace(sys_index_list(ii),sys_index_list(ii+1)-1,Ninterp);
        for id_x = 1 : Nx
            for id_y = 1 : Ny
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

    signal_one_cycle = zeros(size(one_cycle_video, 3),size(one_cycle_video,4),2);
    signal_one_cycle(:,:,1) = squeeze(sum(one_cycle_video .* mask, [1 2])/nnz(mask));

%     figure(33)
%     for ii = 1 : size(one_cycle_video, 4)
%         plot( ...
%             (1:Ninterp),signal_one_cycle(:, ii),'-k', ...
%             'LineWidth',2) ;
%         hold on
%     end
%     title('average blood flow velocity estimate in in-plane retinal arteries');
%     legend(' arterial pulse');
%     fontsize(gca,12,"points") ;
%     ylabel('blood flow velocity (mm/s)');
%     pbaspect([1.618 1 1]) ;
%     set(gca, 'LineWidth', 2);
%     axis tight;
%     hold off

    tmp = zeros(size(video,1), size(video,2), Ninterp);
    ctr = 0;
    for ii = 1:M
%         if (dataReliabilityIndex(ii) < 50)
        if (dataReliabilityIndex(ii) > Reliability_index_threshold)
            tmp = tmp + squeeze(one_cycle_video(:,:,:,ii)); % average selected cycles
            ctr = ctr + 1;
            signal_one_cycle(:,ii,2) = ones(1,size(signal_one_cycle,1));

        end
    end
    % selectedPulseIdx = find(dataReliabilityIndex(ii) < 50);

    if(ctr == 0)
        tmp = squeeze(mean(one_cycle_video,4));% average all cycles
    else
        tmp = tmp /ctr;
    end
else % if M = 0
    1;
end

% FIXME : create zerolag correlation matrix between all zero-mean pulses,
% then svd, to identify the most correlated pulses

% FIXME : export selected pulse indexes

selectedPulseIdx = 0;



one_cycle_video = tmp; % average all detected cycles
oneP = squeeze(sum(one_cycle_video .* mask,[1 2]) / nnz(mask));
[min_val,shift] = min(oneP); % find bottom systole
one_cycle_video = circshift(one_cycle_video, -shift, 3);
signal_one_cycle = circshift(signal_one_cycle,-shift,1); % shift pulse to start & end with bottom diastole

disp('done.');
