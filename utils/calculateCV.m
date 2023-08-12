function [cv,cvs] = calculateCV(spks)

%% Caculates the Coefficient of Variation

% Should have input of nTrials X time

% spikeCounts = sum(spks,2); % Get spike counts across time (dim 2)
FRs = firingRate(spks,'bin_size',100);

cvs = std(FRs,0,2)./mean(FRs,2);

cv = nanmean(cvs);
