function [ARI, ARImap] = construct_resistivity_index(v_RMS, maskArtery)
    % https://en.wikipedia.org/wiki/Arterial_resistivity_index
    % arterial resistivity

    v_RMS = double(v_RMS);

    arterialPulse = squeeze(sum(v_RMS .* maskArtery, [1 2]));
    arterialPulse = arterialPulse / nnz(maskArtery);
    % arterialPulse = smoothdata(arterialPulse, 'rlowess');

    %% avg. arterial resistivity index

    [maxAP, maxAPidx] = max(arterialPulse(:));
    [minAP, minAPidx] = min(arterialPulse(:));

    figure, plot(arterialPulse), hold on
    scatter(maxAPidx, maxAP);
    scatter(minAPidx, minAP);
    yline(maxAP), yline(minAP)
    hold off
    ARI = (maxAP - minAP) / maxAP;

    %% arterial resistivity map values

    ARImap = squeeze((v_RMS(:, :, maxAPidx) - v_RMS(:, :, minAPidx)) ./ v_RMS(:, :, maxAPidx));
    ARImap(ARImap > 1) = 1;
    ARImap = ARImap .* (ARImap .* maskArtery > 0);
    ARImap(isnan(ARImap)) = 0;

end
