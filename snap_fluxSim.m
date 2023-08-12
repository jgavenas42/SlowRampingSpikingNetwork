%% Simulate Networks

clear all

addpath(genpath('./utils'))

%% Initialize parameters
savepath = 'data/flux_eegfix';

N = 400;
Nclusters = 4;
Ksyns = [0,0.5];


noisescale = 0.5;


Ree_params = readtable('REE_params.csv');

rng_list = Ree_params.RNG;
Ree_list = Ree_params.Ree_flux;
% Ree_list = Ree_params.Ree_ignite;

Ntrials = 100;
tsteps = 10000;
                            
%% Initialize parameters in model    


for ii = 1:length(rng_list)
    for jj = 1:length(Ksyns)
        tic
        connMat_rng = rng_list(ii);
        Ree = Ree_list(ii);
        Ksyn = Ksyns(jj);


        % Create the connectivity matrix
        W = connMatrix('N',N,...
            'Nclusters',Nclusters,...
            'Ree',Ree,...
            'Jee',0.2,...
            'Jclust',1.9,...
            'connMat_rng',connMat_rng);

        % Pre-allocate the spike and synaptic potential variables
        spks = zeros(Ntrials,N,tsteps); % Spiking behavior, Ntrials X Nneurons X time
        Rs = zeros(size(spks)); 
        
        tic
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
        
        inh_FRs = squeeze(mean(spks_sm(:,0.8*N:end,:),2));
        inh_FRs = inh_FRs(:,1:50:end)*400;
        
        PSPs = zeros(size(Rs));
        for j=1:Ntrials
            PSPs(j,:,:) = squeeze(W * squeeze(Rs(j,:,:)));
        end
        
        EEG = squeeze(mean(PSPs(:,1:0.8*N,:),2));
        EEG_inh = squeeze(mean(PSPs,2));
        
        
        % Clear out some large variables to save some data.
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

        tmp = cell(Ntrials,1);
        ig_times_all = nan(Ntrials,1);

        for kk=1:Ntrials
            talign_idxs{kk} = strfind(max_Norm_FRs_binary(kk,:),[repmat([0],base_period,1); repmat([1],th_period,1)]');
            if ~isempty(talign_idxs{kk}),ig_times(kk) = talign_idxs{kk}(1);end 

            tmp{kk} = strfind(max_Norm_FRs_binary(kk,:),[repmat([0],5,1); repmat([1],th_period,1)]');
            if ~isempty(tmp{kk}),ig_times_all(kk) = tmp{kk}(1);end 
        end
        
        
        % Extract threshold-aligned data that has lower resolution!
        % Pre-allocation
        
        talign_FRs = nan(Ntrials,Nclusters,base_period+end_period+th_period);
        talign_inh_FRs = nan(Ntrials,base_period+end_period+th_period);
        for k=1:Ntrials
            if ~isnan(ig_times(k))
                if ig_times(k) + base_period + end_period + th_period-1 <= size(cluster_FRs,3)
                    talign_FRs(k,:,:) = cluster_FRs(k,:,ig_times(k):ig_times(k)+base_period+end_period+th_period-1);
                    talign_inh_FRs(k,:) = inh_FRs(k,ig_times(k):ig_times(k)+base_period+end_period+th_period-1);
                    
                end
            end
        end
        
        time_talign_fr = linspace(-50*base_period/1000,50*end_period/1000,size(talign_FRs,3));
        plot(time_talign_fr,squeeze(talign_FRs(:,max_clust,:))'),title('FRs for Max Clust, aligned')
        
        inc_clust_frs = squeeze(nanmean(talign_FRs(:,max_clust,:),1))';
        min_clust = 1:4;
        min_clust = min_clust(min_clust ~= max_clust);
        dec_clust_frs = squeeze(nanmean(talign_FRs(:,min_clust,:),1))';
        
        inh_clust_frs = squeeze(nanmean(talign_inh_FRs,1));
        

        %% Threshold alignment for EEG and single-unit activities
        ig_times_full = ig_times * 50 - 25;
        talign_EEG = nan(Ntrials,(base_period+end_period)*50);
        talign_EEG_inh = nan(Ntrials,(base_period+end_period)*50);
        talign_spks = nan(Ntrials,N,(base_period+end_period)*50);


        for k=1:Ntrials

            if ~isnan(ig_times(k))
                if ig_times_full(k) + (base_period + end_period)*50-1 <= size(EEG,2)
                    talign_EEG(k,:) = EEG(k,ig_times_full(k):ig_times_full(k)+(base_period+end_period)*50-1);
                    talign_EEG_inh(k,:) = EEG_inh(k,ig_times_full(k):ig_times_full(k)+(base_period+end_period)*50-1);

                    talign_spks(k,:,:) = spks(k,:,ig_times_full(k):ig_times_full(k)+(base_period+end_period)*50-1); 
                end
            end


        end
        
        % Baseline the EEG
        talign_EEG = talign_EEG - mean(talign_EEG(:,1:1000),2);
        talign_EEG_inh = talign_EEG_inh - mean(talign_EEG_inh(:,1:1000),2);
        time_talign = linspace(-50*base_period,50*end_period,size(talign_EEG,2));
        
        
        % Create a data variable for saving the values
        data = [];
        
        data.Ree = Ree;
        data.W = W;
        data.tx = ig_times_full + base_period;
        data.Ksyn = Ksyn;
        data.network = connMat_rng;
        
        data.time = time_talign;
        data.time_fr = time_talign_fr;
        
        data.incFRs = inc_clust_frs;
        data.decFRs = dec_clust_frs;
        data.inhFRs = inh_clust_frs;
        
        
        % EEG data
        data.EEG = nanmean(smoothdata(talign_EEG,2,'gaussian',200),1);
        data.EEG_inh = nanmean(smoothdata(talign_EEG_inh,2,'gaussian',200),1);
        
        data.spks = talign_spks(~isnan(ig_times),1:0.8*N,:);
        data.max_clust = max_clust;
        data.min_clust = min_clust;

        
        t_savepath = [savepath '/Ksyn_' num2str(Ksyn)];
        
        if ~exist(t_savepath,'dir')
            mkdir(t_savepath)
        end

        t_fname = sprintf([t_savepath,'/Network_%d.mat'],connMat_rng);

        
        disp(['saving to: ',t_fname])

        save(t_fname,'-struct','data')
        toc
    end
end
        
