function [s,cvs] = calculateSilhouetteScore(spks,k)

%% Caculates the Silhouette score (i.e. how separable 2 spikings are

% Should have input of nTrials X time

% spikeCounts = sum(spks,2); % Get spike counts across time (dim 2)


ss = zeros(size(spks,1),1);

for ii=1:size(spks,1)
    
    ISIs = diff(find(spks(1,:)==1));

    clust = kmeans(ISIs',k);
%     silhouette(ISIs',clust);
    s = silhouette(ISIs',clust);
    
    ss(ii) = nanmean(s,'all');
end

s = nanmean(ss,'all');
