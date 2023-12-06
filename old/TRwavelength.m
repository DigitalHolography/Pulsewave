function TRwavelength(signal)
    % choose a point r0
    final = zeros(size(signal, 1), size(signal, 1));
    pp = zeros(2*size(signal, 1) - 1, size(signal, 2));
    f0 = 150;

    TRsignal = zeros(size(signal));
    for idx = 0 : size(signal, 1) - 1
        TRsignal(idx + 1, :) = signal(end - idx, :);
    end

%     for r0 = 1 : size(signal, 1)
    r0 = 150;
        for frame = 1 : size(signal, 2)
            pp(:, frame) = xcorr(TRsignal(:,frame), signal(:,f0));
        end
%         final(r0, :) = squeeze(mean(pp, 2));
%     end
    signal_final = squeeze(mean(pp, 2));
    [~,~,width,~] = findpeaks(signal_final,1:length(signal_final),'MinPeakProminence',(max(signal_final)-min(signal_final))/2);
    disp(width);
    figure(332)
    plot(signal_final);
    
    

end