function [oneCycleVideo, selectedPulseIdx, cyclesSignal, oneCycleVideoM0] = createOneCycle(video, videoM0, mask, sys_index_list, Ninterp, path, ToolBox)
    %   one_cycle() : identifies pulse cycles and average them to one video
    %   sys_index_list : list of systole indexes in video

    PW_params = Parameters_json(path);

    disp('interpolate, average, and shift');
    disp('arterial pulse waveforms');
    disp('over one cardiac cycle...');

    M = length(sys_index_list) - 1; % M : pulse #
    oneCycleVideo = zeros(size(video, 1), size(video, 2), Ninterp, M);
    oneCycleVideoM0 = zeros(size(video, 1), size(video, 2), Ninterp, M);


    Nx = size(video, 1);
    Ny = size(video, 2);

    if M > 0 % we have detected at least two systoles
        %calculate the signal
        cyclesSignal = zeros(M, Ninterp);
        peak2Peak = zeros(1, M);
        meanSignal = zeros(1, M); 

        for ii = 1:M
            interp_range = linspace(sys_index_list(ii), sys_index_list(ii + 1) - 1, Ninterp);
            frame_range = sys_index_list(ii):sys_index_list(ii + 1) - 1;
            tmp = squeeze(sum(video(:, :, frame_range) .* mask, [1 2]) / nnz(mask));
            %        [idxOutPw, tmp] = discardPulseWaveOutliers( tmp,3);
            %       tmp= smoothdata( tmp,'lowess');
            cyclesSignal(ii, :) = interp1(frame_range, tmp, interp_range);
            peak2Peak(ii) = max(cyclesSignal(ii, :)) - min(cyclesSignal(ii, :));
            % meanSignal(ii) = mean(meanSignal(ii, :));
        end
        
        medianPeak2Peak = median(peak2Peak);
        cyclesSignal = cyclesSignal ./ peak2Peak' .* medianPeak2Peak;

        %check pulse integrity
        pulseMatrix = zeros(M, Ninterp);
        %     cleanSignal = squeeze(mean(one_cycle(:,:,:,:),[1 2 4]));
        cleanSignal = squeeze(mean(cyclesSignal, 1));
        cleanSignal = cleanSignal - mean(cleanSignal(:));

        % FIXME
        %     one_cycle = one_cycle_video(:,:,:,:) .* mask;
        %     one_cycle = one_cycle_video(:,:,:,:);
        dataReliabilityIndex = zeros(M, 1);

        %     figure(34)
        %     plot(cleanSignal, LineWidth=2);

        %     global_noise = zeros(size(cycles_signal));
        %     for ii = 1:M
        %         %         global_noise(ii, :) = sqrt(abs(abs(cycles_signal(ii, :)).^2 - abs(cleanSignal).^2));
        %         global_noise(ii, :) = (abs(cycles_signal(ii, :) - cleanSignal));
        %     end
        %     global_noise_std = std(reshape(global_noise, Ninterp * M, 1));

        for ii = 1:M % for each detected pulse, loop
            %         signal = squeeze(mean(one_cycle(:,:,:,ii),[1 2]));
            signal = cyclesSignal(ii, :);
            signal = signal - mean(signal(:));

            C = sum(signal .* cleanSignal) / Ninterp;
            %
            %         hold on
            %         plot(signal);

            pulseMatrix(ii, :) = signal;
            %     cleanSignal = smoothdata(signal,'lowess');
            %                 noise = sqrt(abs(abs(signal).^2 - abs(cleanSignal).^2));
            %         noise = (abs(signal - cleanSignal));
            %         disp(std(noise))
            %         idxOutNoise = find(noise > 0.5* global_noise_std);
            %         noiseIntensity = squeeze(mean(noise(idxOutNoise)));
            %         idxOutNoise = find(noise>PW_params.oneCycle_outNoiseThreshold*std(noise));
            %         dataReliabilityIndex(ii) = ceil(100*(1-(length(idxOutNoise)/length(noise))));
            %         dataReliabilityIndex2(ii) = (100*(1 - (noiseIntensity)/max(abs(cycles_signal), [], "all")));
            dataReliabilityIndex(ii) = C;
            %         disp(['data reliability index for pulse ', num2str(ii), ' : ', num2str(dataReliabilityIndex(ii)), ' %']);
            %         disp(['data reliability index 2 for pulse ', num2str(ii), ' : ', num2str(dataReliabilityIndex2(ii)), ' %']);
        end

        dataReliabilityIndex = dataReliabilityIndex ./ max(dataReliabilityIndex, [], 'all') * 100;

        for ii = 1:M
            disp(['data reliability index for pulse ', num2str(ii), ' : ', num2str(dataReliabilityIndex(ii)), ' %']);
        end

        % [U,S,V] = svd(pulseMatrix);
        % s1 = S(1,1);
        % S = zeros(size(S));
        % S(1,1) = s1;
        % pulseMatrixFiltered = U * S * V';
        % [row col] = find(S>.3*max(S(:)))

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

        tmp = zeros(size(video, 1), size(video, 2), Ninterp);
        tmpM0 = zeros(size(video, 1), size(video, 2), Ninterp);
        %     ctr = 0;
        cycles_accepted = (dataReliabilityIndex > PW_params.oneCycle_outNoiseThreshold);
        %cycles_accepted = [1 ;1 ;0;0];
        %     for ii = 1:M
        % %         if (dataReliabilityIndex(ii) < 50)
        %         if (dataReliabilityIndex(ii) > PW_params.oneCycle_dataReliabilityThreshold)
        % %             tmp = tmp + squeeze(one_cycle_video(:,:,:,ii)); % average selected cycles
        % %             ctr = ctr + 1;
        % %             signal_one_cycle(:,ii,2) = ones(1,size(signal_one_cycle,1));
        %               cycles_reliability(ii) = 1;
        %         end
        %     end
        % selectedPulseIdx = find(dataReliabilityIndex(ii) < 50);

        single_cycles = zeros(Nx, Ny, Ninterp, M);
        single_cyclesM0 = zeros(Nx, Ny, Ninterp, M);
        selectedPulseIdx = find(cycles_accepted);
        %     parfor ii = 1:M % for each detected pulse, loop
        if ~isempty(selectedPulseIdx) %at least one cycle was accepted

            for jj = 1:length(selectedPulseIdx)
                ii = selectedPulseIdx(jj);
                interp_range = linspace(sys_index_list(ii), sys_index_list(ii + 1) - 1, Ninterp);
                frame_range = sys_index_list(ii):sys_index_list(ii + 1) - 1;

                parfor id_x = 1:Nx

                    for id_y = 1:Ny
                        %                 one_cycle_video(id_x,id_y,:,ii) = interp1((sys_index_list(ii):sys_index_list(ii+1)-1),squeeze(video(id_x, id_y,sys_index_list(ii):sys_index_list(ii+1)-1)),interp_range);
                        single_cycles(id_x, id_y, :, ii) = interp1(frame_range, squeeze(video(id_x, id_y, frame_range)), interp_range);
                        single_cyclesM0(id_x, id_y, :, ii) = interp1(frame_range, squeeze(videoM0(id_x, id_y, frame_range)), interp_range);
                    end

                end

            end

            tmp = squeeze(sum(single_cycles, 4)) / length(selectedPulseIdx);
            tmpM0 = squeeze(sum(single_cyclesM0, 4)) / length(selectedPulseIdx);
        else
            1; % never happens because of the normalization
        end

    else % if M = 0
        fprintf("Less than two systoles were detected");
        cyclesSignal = squeeze(sum(video(:, :, :) .* mask, [1 2]) / nnz(mask));
        cyclesSignal = interp1(cyclesSignal, 1:Ninterp);
        tmp = interp3(video, size(video, 1), size(video, 2), Ninterp);
        tmpM0 = interp3(videoM0, size(video, 1), size(video, 2), Ninterp);
        selectedPulseIdx = 0;

    end

    % FIXME : create zerolag correlation matrix between all zero-mean pulses,
    % then svd, to identify the most correlated pulses
    oneCycleVideo = tmp; % average all detected cycles
    oneCycleVideoM0 = tmpM0; % average all detected cycles

    oneP = squeeze(sum(oneCycleVideo .* mask, [1 2]) / nnz(mask));
    [minVal, shift] = min(oneP); % find bottom systole
    oneCycleVideo = circshift(oneCycleVideo, -shift, 3);
    oneCycleVideoM0 = circshift(oneCycleVideoM0, -shift, 3);
    cyclesSignal = circshift(cyclesSignal, -shift, 2); % shift pulse to start & end with bottom diastole

    disp('done.');