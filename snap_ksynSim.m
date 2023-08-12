%% Simulate Networks

clear all

addpath(genpath('./utils'))

%% Initialize parameters

nettype = 'flux';
savepath = ['data/2clusters/' nettype];

N = 400;
Nclusters = 2;
Ksyns = 0:0.1:0.5;


noisescale = 0.5;


Ree_params = readtable('REE_params.csv');

rng_list = Ree_params.RNG;
Ree_list = Ree_params.Ree_flux * 0.6;
% Ree_list = Ree_params.Ree_ignite;

Ntrials = 100;
tsteps = 10000;
                            
%% Initialize parameters in model    


for ii = 1:length(rng_list)
    connMat_rng = rng_list(ii);
    Ree = Ree_list(ii);
    % Create the connectivity matrix
        W = connMatrix('N',N,...
            'Nclusters',Nclusters,...
            'Ree',Ree,...
            'Jee',0.2,...
            'Jclust',1.9,...
            'connMat_rng',connMat_rng);
    
    for jj = 1:length(Ksyns)
        tic
        Ksyn = Ksyns(jj);

        % Pre-allocate the spike and synaptic potential variables
        spks = zeros(Ntrials,N,tsteps); % Spiking behavior, Ntrials X Nneurons X time
        Rs = zeros(size(spks));
        
        % Run Simulation
        parfor kk=1:Ntrials
            [spks(kk,:,:), Rs(kk,:,:)] = snap_singletrial(W,'tsteps',tsteps,...
                'noisescale',noisescale,...
                'Ksyn',Ksyn); 
        end
        
        disp(['simulation done'])
        toc
        
        spks_sm = smoothdata(spks,3,'gaussian',400);
        cluster_FRs = squeeze(mean(reshape(spks_sm(:,1:0.8*N,:),Ntrials,0.8*N/Nclusters,Nclusters,[]),2));
        cluster_FRs = cluster_FRs(:,:,1:50:end)*400;
        
        PSPs = zeros(size(Rs));
        for j=1:Ntrials
            PSPs(j,:,:) = squeeze(W * squeeze(Rs(j,:,:)));
        end
        
        EEG = squeeze(mean(PSPs(:,1:0.8*N,:),2));
        
        clear Rs PSPs
        
        %% Threshold Alignment
        
        % Pre-processing for threshold alignment
        base_period = 60; % how many continuous 50 ms bins "sub-threshold" (20 = 1 sec)
        end_period = 10; % how many continuous 50 ms bins to retain at the end? (10 = 500 ms)
        th_period = 1; % how many continuous 50 ms bins "above threshold"
        [t_max,max_clust] = max(mean(cluster_FRs,[1 3]));

        % Select only the cluster with maximum firing rate
        max_cluster_FRs = squeeze(cluster_FRs(:,max_clust,:));

        % Normalize the firing rate similar to Fried et al
        max_Norm_FRs = (max_cluster_FRs - mean(max_cluster_FRs(:,1:10),2)) ./ (max(max_cluster_FRs,[],2));
        % Binarize to 0.5 normalized value, also similar to Fried et al
        % findings
        max_Norm_FRs_binary = max_Norm_FRs > 0.5;
        % max_Norm_FRs_binary = max_Norm_FRs > 1.96;
        time_fr = linspace(0,tsteps/1000,size(cluster_FRs,3));
                        
            
        % Loop through trials to find threshold-crossings
        % What we want is a vector that has "baseline" activity for at
        % least 2.5 seconds and then hits some threshold... so create a
        % vector and match an appropriate vector to the binary vector.
        talign_idxs = cell(Ntrials,1);
        ig_times = nan(Ntrials,1);

        for kk=1:Ntrials
            talign_idxs{kk} = strfind(max_Norm_FRs_binary(kk,:),[repmat([0],base_period,1); repmat([1],th_period,1)]');
            if ~isempty(talign_idxs{kk}),ig_times(kk) = talign_idxs{kk}(1);end 
        end
        
        % Now iterate over trials using the ig_times and grab spiking data
        ig_times_full = ig_times * 50 - 25;        
        talign_spks = nan(Ntrials,N,(base_period+end_period)*50);
        time_talign = linspace(-50*base_period,50*end_period,size(talign_spks,3));
        talign_EEG = nan(Ntrials,(base_period+end_period)*50);
        
        for k=1:Ntrials
            if ~isnan(ig_times(k))
                if ig_times_full(k) + (base_period + end_period)*50-1 <= size(spks,3)
                    talign_EEG(k,:) = EEG(k,ig_times_full(k):ig_times_full(k)+(base_period+end_period)*50-1);
                    talign_spks(k,:,:) = spks(k,:,ig_times_full(k):ig_times_full(k)+(base_period+end_period)*50-1); 
                end
            end
        end
        
        talign_EEG = talign_EEG - mean(talign_EEG(:,1:1000),2);
        
        %% Threshold-alignment of single-neuron data
        talign_spks = talign_spks(~isnan(ig_times),1:0.8*N,:);

        
        
        type = cell(size(talign_spks,2),1);
    
        if size(talign_spks,1) < 10 % If less than 10 trials...
            type(:) = {'bad'};
            incSpks = nan(3500,1);
            decSpks = nan(3500,1);
            nonSpks = nan(3500,1);
            acfs_all = nan(size(spks,2),1001);
            acfs_clust = nan(4,1001);
        else
            acfs_all = zeros(size(spks,2),1001);

            for kk = 1:size(talign_spks,2)
                spks_this = squeeze(talign_spks(:,kk,:));

                idxs = ((time_talign >= -3000) & (time_talign < -2000));
                baseFRs = mean(spks_this(:,idxs),2)*1000;

                idxs = ((time_talign >= -400) & (time_talign < 0));
                stopFRs = nanmean(spks_this(:,idxs),2)*1000;

                % Rank-sum test to classify
                [P,H,STATS] = ranksum(stopFRs,baseFRs,'alpha',0.01);

                if isnan(H) || ~isfield(STATS,'zval')
                    type{kk} = 'bad';
                elseif ~H
                    type{kk} = 'nochange';
                elseif (STATS.zval > 0)
                    type{kk} = 'inc';
                elseif (STATS.zval < 0)
                    type{kk} = 'dec';
                end
                
                spks_this_sm = smoothdata(spks_this,2,'movmean',10);
                spks_this_sm = spks_this_sm(~isnan(spks_this_sm(:,1)),:);

                acfs_this = zeros(size(spks_this_sm,1),1001);
                for ll = 1:size(spks_this_sm,1)
                    acfs_this(ll,:) = autocorr(spks_this_sm(ll,:)',1000);
                end
                acfs_all(kk,:) = nanmean(acfs_this,1);
            end
            spks_sm = smoothdata(talign_spks,3,'gaussian',400)*1000;

            nInc = sum(strcmp(type,'inc'));
            nDec = sum(strcmp(type,'dec'));
            nNon = sum(strcmp(type,'non'));

            incSpks = squeeze(nanmean(spks_sm(:,strcmp(type,'inc'),:),[1 2]));
            decSpks = squeeze(nanmean(spks_sm(:,strcmp(type,'dec'),:),[1 2]));
            
            
            spks_sm = smoothdata(talign_spks,3,'movmean',10);
            cluster_FRs = squeeze(mean(reshape(spks_sm,size(talign_spks,1),80,4,[]),2));

            cluster_FRs = cluster_FRs(~isnan(cluster_FRs(:,1,1)),:,:);
            acfs_clust = zeros(4,1001);
            for ll = 1:size(cluster_FRs,2)
                frs_this = squeeze(cluster_FRs(:,ll,:));
                acfs_this = zeros(size(frs_this,1),1001);
                for kk = 1:size(frs_this,1)
                    acfs_this(kk,:) = autocorr(frs_this(kk,:)',1000);
                end
                acfs_clust(ll,:) = nanmean(acfs_this,1);
            end
            
        end
        
        
        
        
        % Create a data variable for saving the values
        data = [];
        
        data.Ree = Ree;
        data.Ksyn = Ksyn;
        data.network = connMat_rng;
        
        data.time = time_talign;
        
        data.incSpks = incSpks;
        data.decSpks = decSpks;
        
        data.EEG = talign_EEG;
        
        data.acfs_all = acfs_all;
        data.acfs_clust = acfs_clust;
        
        
        
        
        if ~exist(savepath,'dir')
            mkdir(savepath)
        end

        t_fname = sprintf([savepath,'/Ksyn_%d_Network_%d.mat'],Ksyn,connMat_rng);

        
        disp(['saving to: ',t_fname])

        save(t_fname,'-struct','data')
        toc
    end
end
        
