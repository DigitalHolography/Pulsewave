function [interpSig,avgLength,interpstdSig] = interpSignal(sig,findIndxs,Ninterp,stdSig)
    % Retourne une version interpolée du signal sur plusieurs cycles
    % désignés par findIndxs
    signal = sig(findIndxs(1):findIndxs(end));
    avgLength = mean(findIndxs(2:end)-findIndxs(1:end-1));
    n = length(signal);


    interpSig = zeros([1 Ninterp]);
    for i=1:(length(findIndxs)-1)
        X=findIndxs(i):findIndxs(i+1);
        interpSig = interpSig + interp1(X,sig(X),linspace(X(1),X(end),Ninterp));
    end
    interpSig = interpSig/(length(findIndxs)-1);
    
    if nargin > 3
        interpstdSig = zeros([1 Ninterp]);
        for i=1:(length(findIndxs)-1)
            X=findIndxs(i):findIndxs(i+1);
            interpstdSig = interpstdSig + interp1(X,stdSig(X),linspace(X(1),X(end),Ninterp)).^2;
        end
        interpstdSig = sqrt(interpstdSig/(length(findIndxs)-1));
    end
end