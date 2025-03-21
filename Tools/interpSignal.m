function [interpSig, avgLength, interpstdSig] = interpSignal(sig, findIndxs, Ninterp, stdSig)
% Retourne une version interpolée du signal sur plusieurs cycles
% désignés par findIndxs

interpSig = zeros([1 Ninterp]);

if length(findIndxs) < 2 %(less than one systole) interpolate all sig
    X = 1:length(sig);
    interpSig = interp1(X, sig(X), linspace(X(1), X(end), Ninterp));
    avgLength = length(sig);

    if nargin > 3
        interpstdSig = zeros([1 Ninterp]);
        X = 1:length(sig);
        interpstdSig = interpstdSig + interp1(X, stdSig(X), linspace(X(1), X(end), Ninterp)) .^ 2;
    end

    return
end

avgLength = mean(findIndxs(2:end) - findIndxs(1:end - 1));

for i = 1:(length(findIndxs) - 1)
    X = findIndxs(i):findIndxs(i + 1);
    interpSig = interpSig + interp1(X, sig(X), linspace(X(1), X(end), Ninterp));
end

interpSig = interpSig / (length(findIndxs) - 1);

if nargin > 3
    interpstdSig = zeros([1 Ninterp]);

    for i = 1:(length(findIndxs) - 1)
        X = findIndxs(i):findIndxs(i + 1);
        interpstdSig = interpstdSig + interp1(X, stdSig(X), linspace(X(1), X(end), Ninterp)) .^ 2;
    end

    interpstdSig = sqrt(interpstdSig / (length(findIndxs) - 1));
end

end
