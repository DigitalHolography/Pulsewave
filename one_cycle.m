function one_cycle_video = one_cycle(video, Ninterp)
    videoObject = VideoReader(video);
    h=videoObject.height;
    w=videoObject.width;
    
    NumberFrames=videoObject.NumFrames;
    
    %Video retrieval frame by frame
    I = vid2frames(videoObject,NumberFrames,h,w);
    for pp = 1:size(I, 3)
        I(:,:,pp) = I(:,:,pp) ./ mean(I(:,:,pp), [1 2]);
    end
    mask = std(I, 0, 3);
    mask = imbinarize(im2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
    figure(42)
    imagesc(mask);

%     pulse = zeros(1, 2 * size(I, 3) - 1);
%     for pp = 1:size(I, 3)-1
%         I(:,:,pp) = I(:,:,pp) .* mask;
%         I(:,:,pp+1) = I(:,:,pp+1) .* mask;
%         J(:,:) = (I(:,:,pp) + I(:,:,pp + 1)) ./ 2;
%         pulse(2*pp - 1) = squeeze(mean(I(:,:,pp), [1 2]));
%         pulse(2*pp) = squeeze(mean(J(:,:), [1 2]));
%     end
%     pulse(2 * size(I, 3) - 1) = squeeze(mean(I(:,:,size(I, 3)), [1 2]));
    pulse = squeeze(mean(I .* mask, [1 2]));
    figure(1);
    plot(pulse);


    pulse_init = pulse - mean(pulse, "all");
    C = I;
    for kk = 1:size(I, 3)
        C(:,:,kk) = I(:,:,kk) - squeeze(mean(I, 3));
    end
    pulse_init_3d = zeros(size(I));
    for mm = 1:size(I, 1)
        for pp = 1:size(I, 2)
            pulse_init_3d(mm,pp,:) = pulse_init;
        end
    end
    C = C .* pulse_init_3d;
    C = squeeze(mean(C, 3));
    figure(43)
    imagesc(C);
    colorbar

    mask = C > max(C(:))*0.2;
    figure(2)
    imagesc(mask);
    pulse = squeeze(mean(I .* mask, [1 2]));
    figure(1);
    plot(pulse);


    pulse_init = pulse - mean(pulse, "all");
%     C = I - mean(I, 3);
%     pulse_init_3d = zeros(size(I));
%     for mm = 1:size(I, 1)
%         for pp = 1:size(I, 2)
%             pulse_init_3d(mm,pp,:) = pulse_init;
%         end
%     end
%     C = C .* pulse_init_3d;
%     C = squeeze(mean(C, 3));
%     figure(43)
%     imagesc(C);
%     colorbar
% 
%     tic
%     I0 = I - mean(I, 3);
%     for kk = 1:size(I, 3)
%         R = I0 .* circshift(pulse_init_3d, kk, 3);
%         C(:,:,kk) = squeeze(mean(R, 3));
%     end
%     C = imgaussfilt3(C, 2);
%     [max_C_3, id_max] = max(C, [], 3);
%     figure(45);
%     imagesc(id_max);
%     toc


%     derivative = diff(pulse)./diff(1:length(pulse));
%     derivative_mask = derivative > mean(derivative);
%     systole_begin_index = [];
%     for ii = 2:length(derivative_mask)
%         if (and(derivative_mask(ii), ~derivative_mask(ii - 1)))
%             systole_begin_index(end + 1) = ii;
%         end
%     end
%     figure(2);
%     hold on
%     for ii = 2:length(systole_begin_index)
%         plot(pulse(systole_begin_index(ii-1):systole_begin_index(ii)));
%     end
%     hold off
    pulse_init = detrend(pulse_init);
    Adiff = diff(pulse_init);
    x=1:length(Adiff);
    [pks,lsor] = findpeaks(Adiff,x,'MinPeakDistance',35,'SortStr','descend');
    Npeaks = length(lsor);
    lsor = sort(lsor(1:Npeaks));
    M = Npeaks-1;
    A1 = zeros(Ninterp,M);
    for ii = 1:M
        interp_range = linspace(lsor(ii),lsor(ii+1)-1,Ninterp); %range of interpolation for each cycle
        A1(:,ii) = interp1((lsor(ii):lsor(ii+1)-1),pulse_init(lsor(ii):lsor(ii+1)-1),interp_range); %interpolation of each cycle
    end %ii
    A1 = circshift(A1,30,1); %shifts the beginning of the cycle at the beginning og the matrice

    A1_video = zeros(size(I,1), size(I,2), Ninterp, M);

    for mm = 1:M
        for id_x = 1 : size(I,1)
            for id_y = 1 : size(I,2)
                interp_range = linspace(lsor(ii),lsor(ii+1)-1,Ninterp);
                A1_video(id_x,id_y,:,mm) = interp1((lsor(ii):lsor(ii+1)-1),squeeze(I(id_x, id_y,lsor(ii):lsor(ii+1)-1)),interp_range);
            end
        end
    end

    one_cycle_video = mat2gray(mean(A1_video, 4));

%     w = VideoWriter('one_cycle_video.avi');
%     open(w);
%     for i = 1:size(one_cycle_video,3)
%         writeVideo(w, one_cycle_video(:,:,i));   
%     end
%     close(w)
end