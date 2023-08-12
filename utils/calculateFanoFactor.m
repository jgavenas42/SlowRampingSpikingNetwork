function [ff,ffs] = calculateFanoFactor(spks)

%% Caculates the Fano Factor

% Should have input of nTrials X time

% spikeCounts = sum(spks,2); % Get spike counts across time (dim 2)
ffs = zeros(size(spks,1),1);

for ii=1:size(spks,1)
    
    ISIs = diff(find(spks(ii,:)==1));

    ffs(ii) = var(ISIs)./(mean(ISIs)^2);
end

ff = nanmean(ffs);


