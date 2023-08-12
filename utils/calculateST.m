function ST = calculateST(spks,Nclusters)

%% Quantifies the presence of slow-switching assemblies
% As in Schaub et al., 2015, we calculate the spiking variability over time
% in the cluster of firing rates according to spikes.


FRs = firingRate(spks,'bin_size',100);

cluster_FRs = squeeze(mean(reshape(FRs,size(FRs,1)/Nclusters,Nclusters,[]),1));

ST_real = mean(std(cluster_FRs,0,2));

ST_shuffle = zeros(20,1);
for ii=1:20
    FRs_shuffle = FRs(randperm(size(FRs, 1)), :);
    cluster_FRs = squeeze(mean(reshape(FRs_shuffle,size(FRs_shuffle,1)/Nclusters,Nclusters,[]),1));
    ST_shuffle(ii) = mean(std(cluster_FRs,0,2));
end

ST = ST_real - mean(ST_shuffle);
