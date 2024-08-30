function [API, APImap] = construct_pulsatility_index(v_RMS, maskArtery)

    % arterial resistivity

    v_RMS = double(v_RMS);

    arterialPulse = squeeze(sum(v_RMS .* maskArtery, [1 2]));
    arterialPulse = arterialPulse / nnz(maskArtery);

    %% avg. arterial resistivity index

    [maxAP, maxAPidx] = max(arterialPulse(:));
    [minAP, minAPidx] = min(arterialPulse(:));
    meanAP = mean(arterialPulse);

    API = (maxAP - minAP) / meanAP;

    %% arterial resistivity map values

    APImap = squeeze((v_RMS(:, :, maxAPidx) - v_RMS(:, :, minAPidx)) ./ mean(v_RMS, 3));
    APImap(APImap > 1) = 1;
    APImap = APImap .* (APImap .* maskArtery > 0);
    APImap(isnan(APImap)) = 0;

end
