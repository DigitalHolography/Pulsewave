function SegmentationAV(one_pulse_video, directory, filename)

%     videoObject = VideoReader(one_pulse_video);
%     h=videoObject.height;
%     w=videoObject.width;
%     
%     NumberFrames=videoObject.NumFrames;
%     
%     I = zeros(h,w, NumberFrames);
% 
%     %Video retrieval frame by frame
%     for i = 1 : NumberFrames
%         Frame = read(videoObject, i);
%         I(:,:,i)=rgb2gray(Frame);
%     end
    
    I = one_pulse_video ; 
    % normalize through laser momentary intensity FIXME
    for pp = 1:size(I, 3)
        I(:,:,pp) = I(:,:,pp) ./ mean(I(:,:,pp), [1 2]);
    end
    
    % identify dominant arteries pulse
    mask_arteries = std(I, 0, 3);
    mask_arteries = imbinarize(im2gray(mask_arteries), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.05);
    I_arteries = I .* mask_arteries;
    % calculate the vector of pulse in arteries
    pulse_arteries = squeeze(mean(I_arteries, [1 2]));

%     pulse_init = pulse - mean(pulse, "all");
%     C = I - mean(I, 3);
    C = I;
    pulse_arteries_3d = zeros(size(I_arteries));
    for mm = 1:size(I, 1)
        for pp = 1:size(I, 2)
            pulse_arteries_3d(mm,pp,:) = pulse_arteries;
        end
    end

    tic
    I0 = I - mean(I, 3);
    for kk = 1:size(I, 3)
        R = I0 .* circshift(pulse_arteries_3d, kk, 3);
        C(:,:,kk) = squeeze(mean(R, 3));
    end
%     C = imgaussfilt3(C, 2);
    [max_C_3, id_max] = max(C, [], 3);
    figure(1);
    imagesc(id_max);
    colorbar ; 
    toc
    
    print('-f1','-dpng',fullfile(directory,strcat(filename,'_AVmap.png'))) ; 
    

    % pulses 
    nb_frame = size(one_pulse_video,3) ;
    arteries = zeros(size(id_max)); 
    arteries (id_max >0.95*nb_frame) = 1 ; 
    figure(10)
    imagesc(arteries) ;
    veins = zeros(size(id_max)) ; 
    veins(and(id_max>0.15*nb_frame,id_max<0.3*nb_frame)) = 1 ;  
    figure(11)
    imagesc(veins) ; 

    I_arteries = I .*arteries;
    % calculate the vector of pulse in arteries
    pulse_arteries = squeeze(sum(I_arteries, [1 2]))/nnz(I_arteries);
    pulse_arteries = pulse_arteries - min(pulse_arteries) ; 

    I_veins = I .*veins ; 
    pulse_veins = squeeze(sum(I_veins, [1 2]))/nnz(I_veins);
    pulse_veins = pulse_veins - min(pulse_veins) ; 
    figure(2)
    plot(pulse_arteries,'k') ; 
    hold on 
    plot(pulse_veins,'b') ; 
    hold off
    legend('Arteries','Veins') ; 
    
    print('-f2','-dpng',fullfile(directory,strcat(filename,'_pulse_veins.png'))) ; 


end