%% SNAP talign

% Spare work for the threshold alignment functions...
% Load a function where things were simulated
clear all
close all

addpath(genpath('./utils'))


% fpath = './data_th/Ksyn_0.75/';
fpath = '/Volumes/GDrive/snap/data_th/Ksyn_0';

files = dir([fpath,'*.mat']);
t_fid = fopen([fpath,'filesplotted2.txt'],'a+');
t_A = fileread([fpath,'filesplotted2.txt']);

for i=1:length(files)
    t_fname = files(i).name;
    t_idx = strfind(t_A,t_fname);
    
    temp = split(files(i).name,'.mat');

    fname = temp{1};
    searchstring = 'Ree';
    
    if contains(t_fname,searchstring) && isempty(t_idx)
        
        if ~exist([fpath,fname],'dir'),mkdir([fpath,fname]),end
        
        load([fpath,t_fname])
        %% Visualizing individual trials
%         start_idx = 5;
%         num_trials = 4;
%         figure('Position', [10 10 900 600])
% 
%         for j=start_idx:start_idx+num_trials
%             subplot(num_trials,1,mod(j,num_trials)+1)
%             plot((1:50:tsteps)/(1000),squeeze(cluster_FRs(j,:,:))'), hold on
%         %     legend;
%             plot((1:50:tsteps)/(1000),mean(squeeze(cluster_FRs(j,:,:)),1),'k','LineWidth',2)
%             xlabel('time (s)'),ylabel('Firing Rate (Hz)')
%             title('FRs on each cluster + avg')
% 
%         %             subplot(i,2,2)
%         %             plot((1:ex_tsteps)/1000,movmean(cluster_PSPs',50)), hold on
%         %             plot((1:ex_tsteps)/1000,movmean(mean(cluster_PSPs,1),50),'k','LineWidth',2)
%         %             xlabel('time (s)'),ylabel('mV (post synaptic potentials)'),title('PSPs on each cluster + avg')
%         end
%         sgtitle([num2str(num_trials),'individual trials of activity'])
% 
%         saveas(gcf,[fpath,fname,'/indivtrials'],'epsc')

        %% Pre-processing for threshold alignment
        base_period = 50; % how many continuous 50 ms bins "sub-threshold" (20 = 1 sec)
        end_period = 10; % how many continuous 50 ms bins to retain at the end? (10 = 500 ms)
        th_period = 1; % how many continuous 50 ms bins "above threshold"
        [t_max,max_clust] = max(mean(cluster_FRs,[1 3]));

        % Select only the cluster with maximum firing rate
        max_cluster_FRs = squeeze(cluster_FRs(:,max_clust,:));

        % Normalize the firing rate similar to Fried et al
        max_Norm_FRs = (max_cluster_FRs - mean(max_cluster_FRs(:,1:10),2)) ./ (max(max_cluster_FRs,[],2));
        % max_Norm_FRs = (max_cluster_FRs - mean(max_cluster_FRs(:,1:10),2)) ./ std(max_cluster_FRs(:,1:10),0,2);
        % Binarize to 0.5 normalized value, also similar to Fried et al
        % findings
        max_Norm_FRs_binary = max_Norm_FRs > 0.5;
        % max_Norm_FRs_binary = max_Norm_FRs > 1.96;
        time_fr = linspace(0,tsteps/1000,size(cluster_FRs,3));
%         figure
%         plot(time_fr,max_Norm_FRs');title('Normalized FRs for Max Clust unaligned')
%         yline(0.5,'--','LineWidth',2)
%         saveas(gcf,[fpath,fname,'/MaxClustDDM'],'epsc')

        %% Threshold-Alignment Loop
        % What we want is a vector that has "baseline" activity for at
        % least 2 seconds and then hits some threshold... so create a
        % vector and match an appropriate vector to the binary vector.
        talign_idxs = cell(Ntrials,1);
        ig_times = nan(Ntrials,1);

        tmp = cell(Ntrials,1);
        ig_times_all = nan(Ntrials,1);

        for k=1:Ntrials
            talign_idxs{k} = strfind(max_Norm_FRs_binary(k,:),[repmat([0],base_period,1); repmat([1],th_period,1)]');
            if ~isempty(talign_idxs{k}),ig_times(k) = talign_idxs{k}(1);end 

            tmp{k} = strfind(max_Norm_FRs_binary(k,:),[repmat([0],5,1); repmat([1],th_period,1)]');
            if ~isempty(tmp{k}),ig_times_all(k) = tmp{k}(1);end 
        end

%         figure
%         histogram(ig_times_all*50/1000);title('Times of First Threshold-Crossing')
%         saveas(gcf,[fpath,fname,'/TXdistribution'],'epsc')


        % num_igs = sum(~isnan(ig_times));
        talign_FRs = nan(Ntrials,Nclusters,51);
        for k=1:Ntrials

            if ~isnan(ig_times(k))
                talign_FRs(k,:,:) = cluster_FRs(k,:,ig_times(k):ig_times(k)+50);
            end


        end

        %% Plotting for Threshold-Aligned Activities
        time_talign_fr = linspace(-2.5,0,size(talign_FRs,3));
        % plot(time_talign_fr,squeeze(talign_FRs(:,max_clust,:))'),title('FRs for Max Clust, aligned')
%         figure
%         for k=1:Nclusters
%             subplot(2,2,k)
%             if k == max_clust
%                 plot(time_talign_fr,squeeze(talign_FRs(:,k,:))','r')
%             else
%                 plot(time_talign_fr,squeeze(talign_FRs(:,k,:))','b')
%             end
%             hold on
%             plot(time_talign_fr,nanmean(squeeze(talign_FRs(:,k,:)),1),'k','linewidth',2)
%             xlabel('Time rel. threshold crossing')
%             ylabel('Cluster-Averaged Firing Rate')
% 
%         end
%         sgtitle(['TAligned Activities for the 4 clusters, ',num2str(sum(~isnan(ig_times))),'/100 trials retained'])
%         saveas(gcf,[fpath,fname,'/TAligned4Clusters'],'epsc')

        %% Post-Synaptic Potentials
        PSPs = zeros(size(Rs));
        for j=1:Ntrials
            PSPs(j,:,:) = squeeze(W * squeeze(Rs(j,:,:)));
        end
        cluster_PSPs = squeeze(mean(reshape(PSPs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));


        %% Threshold-alignment for PSPs
        % the PSPs are larger than the firing rates, because the firing rates have
        % been binned in increments of 50 ms. The PSPs (and spikes, which we'll use
        % later) have not. We want to select the right time for 'ignitions' in the
        % more fine-grained time. So if the 'ignition time' is index 50 for the 

        ig_times_full = ig_times * 50 - 25;
        talign_PSPs = nan(Ntrials,2501);


        for k=1:Ntrials

            if ~isnan(ig_times(k))
                talign_PSPs(k,:) = nanmean(cluster_PSPs(k,:,ig_times_full(k):ig_times_full(k)+2500),2);
            end


        end

        %Baseline
        talign_PSPs = talign_PSPs - mean(talign_PSPs(:,1:1000),2);
        time_talign = linspace(-2.5,0,size(talign_PSPs,2));
%         figure
%         plot(time_talign,talign_PSPs','g')
%         hold on
%         plot(time_talign,nanmean(talign_PSPs,1),'k','linewidth',2)
%         xlabel('Time rel. threshold crossing')
%         ylabel('mV')
%         title('Post-synaptic potentials on excitatory neurons')
%         saveas(gcf,[fpath,fname,'/TAlignedPSPs_all'],'epsc')
% 
%         figure
% 
%         plot(time_talign,nanmean(talign_PSPs,1),'g','linewidth',2)
%         xlabel('Time rel. threshold crossing')
%         ylabel('mV')
%         title('Post-synaptic potentials on excitatory neurons')
%         saveas(gcf,[fpath,fname,'/TAlignedPSPs'],'epsc')

        
        if Ksyn ~= 1
            figure

            plot(time_talign,smoothdata(nanmean(talign_PSPs,1),'movmean',200),'g','linewidth',2)
            xlabel('Time rel. threshold crossing')
            ylabel('mV')
            title('Post-synaptic potentials on excitatory neurons')
            saveas(gcf,[fpath,fname,'/TAlignedPSPs_smooth'],'epsc')
        end


%         %% Replication of FMK method on simulated data
% 
% 
%         data_all = cell(N,1);
% 
%         for j=1:N
% 
% 
% 
%             data = [];
% 
%             spks_full = full(spks{j});
% 
%             TAlignedSpikes = nan(Ntrials,2501);
% 
%             for k = 1:Ntrials
%                 if ~isnan(ig_times(k))
%                     TAlignedSpikes(k,:) = spks_full(k,ig_times_full(k):ig_times_full(k)+2500);
%                 end
%             end
% 
%             data.TAlignedSpikes = TAlignedSpikes;
% 
% 
%             data.cluster=k;
% 
%             data.ntrials = Ntrials;
% 
%             data.time = linspace(-2.5,0,size(TAlignedSpikes,2));
% 
%             idxs = ((data.time >= -2.5) & (data.time < -1.5));
%             baseFRs = mean(TAlignedSpikes(:,idxs),2)*1000;
%         %         baseFRs = mean(BaselineAlignedSpikes,2)*1000;
%             data.baseFRs = baseFRs;
% 
%             stopFRs = nanmean(TAlignedSpikes(:,2001:end),2)*1000;
%             data.stopFRs = stopFRs;
% 
% 
%             % Rank-sum test to classify
%             [P,H,STATS] = ranksum(stopFRs,baseFRs,'alpha',0.01);
% 
%             data.p = P;
%             data.zval = STATS.zval;
%             if j > 0.8*N
%                 data.type = 'inh';
%             elseif isnan(H)
%                 data.type = 'bad';
%             elseif ~H
%                 data.type = 'nochange';
%             elseif (STATS.zval > 0)
%                 data.type = 'inc';
% 
%             elseif (STATS.zval < 0)
%                 data.type = 'dec';
% 
%             end
% 
%             % Sigmoid fitting
%             % Fitting Sigmoid
%             fB = num2str(nanmean(baseFRs));
%             fW = num2str(nanmean(stopFRs));
% 
%             if strcmp(data.type,'inc')
%                 model = strcat('(',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
% %                 model2 = strcat(fB,' + exp((x - t0)/gain)');
% 
%                 fo = fitoptions('Method','NonlinearLeastSquares',...
%                                'Algorithm','Levenberg-Marquardt');
%                 ft = fittype(model,'options',fo);
% 
%                 frs = nanmean(firingRate(TAlignedSpikes(~isnan(ig_times),:)),1);
%                 time_fr = linspace(-3.5,0,length(frs));
%                 idxs = find((time_fr > -2.5));
%                 [fitted,gof,output] = fit(time_fr(idxs)',smoothdata(frs(idxs),'gaussian',4)',ft);
%                 data.sigmoid = fitted;
%                 data.sigmoid_gof = gof;
%                 data.gain = fitted.gain;
%                 plot(fitted,time_fr(idxs),smoothdata(frs(idxs),'gaussian',4))
%                 title(['Type: ',data.type,', P = ',num2str(data.p)])
%                 yline(str2double(fB));
%                 yline(str2double(fW));
%                 xlim([-2.5 0])
%                 legend
%             elseif strcmp(data.type,'dec')
%                 model = strcat('(-',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
%                 model2 = strcat(fB,' - exp((x - t0)/gain)');
%                 fo = fitoptions('Method','NonlinearLeastSquares',...
%                                'Algorithm','Levenberg-Marquardt');
%                 ft = fittype(model,'options',fo);
% 
%                 frs = nanmean(firingRate(TAlignedSpikes(~isnan(ig_times),:)),1);
%                 time_fr = linspace(-3.5,0,length(frs));
%                 idxs = find((time_fr > -2.5));
%                 [fitted,gof,output] = fit(time_fr(idxs)',smoothdata(frs(idxs),'gaussian',4)',ft);
%                 data.sigmoid = fitted;
%                 data.sigmoid_gof = gof;
%                 data.gain = fitted.gain;
%                 plot(fitted,time_fr(idxs),smoothdata(frs(idxs),'gaussian',4))
%                 title(['Type: ',data.type,', P = ',num2str(data.p)])
%                 yline(str2double(fB));
%                 yline(str2double(fW));
%                 xlim([-2.5 0])
%                 legend
% 
%             end
% 
% 
% 
% 
%             data.FRs = firingRate(TAlignedSpikes);
%             t_normFRs = (data.FRs - baseFRs);
% 
%             data.normFRs = t_normFRs ./ max(abs(t_normFRs'))';
%             
% 
%             data.tau = calculateTimescale('spk_mat',TAlignedSpikes(~isnan(ig_times),1:1000));
% 
%         %     subplot(2,1,2)
%         %     plot(time_fr,nanmean(data.normFRs,1))
%             data_all{j} = data;
%         %     pause(2)
%         end
%         % saveas(gcf,[fpath,fname,'/TAlignedPSPs'],'.epsc')
%         save([fpath,fname],'data_all')
% 
% 
%         %% Part Two
%         nInc = 0;
%         nDec = 0;
%         nNon = 0;
%         nInh = 0;
% 
%         for k=1:N
% 
%             if strcmp(data_all{k}.type,'inc')
%                 nInc = nInc + 1;
%             elseif strcmp(data_all{k}.type,'dec')
%                 nDec = nDec + 1;
%             elseif strcmp(data_all{k}.type,'nochange')
%                 nNon = nNon + 1;
%             elseif strcmp(data_all{k}.type,'inh')
%                 nInh = nInh + 1;
%             end
%         end
% 
% 
%         decFRs = zeros(nDec,size(t_normFRs,2));
%         incFRs = zeros(nInc,size(t_normFRs,2));
%         nonFRs = zeros(nNon,size(t_normFRs,2));
%         inhFRs = zeros(nInh,size(t_normFRs,2));
% 
% 
%         decNormFRs = zeros(nDec,size(t_normFRs,2));
%         incNormFRs = zeros(nInc,size(t_normFRs,2));
%         nonNormFRs = zeros(nNon,size(t_normFRs,2));
%         inhNormFRs = zeros(nInh,size(t_normFRs,2));
% 
%         decNormFRs_var = zeros(nDec,size(t_normFRs,2));
%         incNormFRs_var = zeros(nInh,size(t_normFRs,2));
%         nonNormFRs_var = zeros(nNon,size(t_normFRs,2));
%         inhNormFRs_var = zeros(nInh,size(t_normFRs,2));
% 
% 
%         decTaus = zeros(nDec,1);
%         incTaus = zeros(nInc,1);
%         nonTaus = zeros(nNon,1);
%         inhTaus = zeros(nInh,1);
% 
% 
%         decIter = 0;
%         incIter = 0;
%         nonIter = 0;
%         inhIter = 0;
% 
%         incGain = zeros(nInc,1);
%         decGain = zeros(nDec,1);
% 
%         incR2 = zeros(nInc,1);
%         decR2 = zeros(nDec,1);
% 
%         incP = zeros(nInc,1);
%         decP = zeros(nDec,1);
% 
%         % incFitted = cell{nInc,2}; % First column fit objects, second gof object
%         % decFitted = cell{nDec,2}; % First column fit objects, second column gof object
% 
% 
% 
%         for l = 1:N
%         %     load([fpath,filenames(i).name],'normFRs','type','BaselineAlignedSpikes','p')
% 
% 
%             if strcmp(data_all{l}.type,'bad'),continue,end
% 
%             tmpFR = nanmean(data_all{l}.normFRs);
% 
%             tmpFRvar = nanvar(data_all{l}.normFRs);
%             tmpFRfull = nanmean(data_all{l}.FRs);
% 
% 
%             tmpTau = data_all{l}.tau;
% 
%             if strcmp(data_all{l}.type,'inc')
%                 incIter = incIter + 1;
%                 incFRs(incIter,:) = tmpFRfull;
% 
% 
%                 incNormFRs(incIter,:) = tmpFR;
%                 incNormFRs_var(incIter,:) = tmpFRvar;
%                 incTaus(incIter) = tmpTau;
% 
%                 incP(incIter) = data_all{l}.p;
% 
%         %         load([fpath,filenames(i).name],'gain','sigmoid_gof')
%                 incGain(incIter) = data_all{l}.gain;
%                 incR2(incIter) = data_all{l}.sigmoid_gof.rsquare;
% 
% 
% 
%             elseif strcmp(data_all{l}.type,'dec')
%                 decIter = decIter + 1;
% 
%                 decFRs(decIter,:) = tmpFRfull;
% 
% 
%                 decNormFRs(decIter,:) = tmpFR;
%                 decNormFRs_var(decIter,:) = tmpFRvar;
%                 decTaus(decIter,:) = tmpTau;
% 
%                 decP(decIter) = data_all{l}.p;
% 
%         %         load([fpath,filenames(i).name],'gain','sigmoid_gof')
%                 decGain(decIter) = data_all{l}.gain;
%                 decR2(decIter) = data_all{l}.sigmoid_gof.rsquare;
% 
% 
%             elseif strcmp(data_all{l}.type,'nochange')
%                 nonIter = nonIter + 1;
% 
%                 nonFRs(nonIter,:) = tmpFRfull;
% 
% 
%                 nonNormFRs(nonIter,:) = tmpFR;
%                 nonNormFRs_var(nonIter,:) = tmpFRvar;
%                 nonTaus(nonIter,:) = tmpTau;
% 
% 
%             elseif strcmp(data_all{l}.type,'inh')
%                 inhIter = inhIter + 1;
% 
%                 inhFRs(inhIter,:) = tmpFRfull;
% 
% 
%                 inhNormFRs(inhIter,:) = tmpFR;
%                 inhNormFRs_var(inhIter,:) = tmpFRvar;
%                 inhTaus(inhIter,:) = tmpTau;
% 
%             end
%         end
% 
% 
%         %% Plotting
%         % time_fr = -3.45:0.05:0;
%         figure
% 
%         subplot(1,2,1)
%         plot(time_fr,nanmean(incNormFRs),'r')
%         hold on, grid on
%         plot(time_fr,nanmean(decNormFRs),'b')
%         plot(time_fr,nanmean(nonNormFRs),'k')
%         % plot(time_fr,mean(inhNormFRs),'g')
%         title('Normalized Firing Rates')
%         xlim([-2,0])
%         legend('Inc','Dec','Non','location','northwest')
% 
% 
%         subplot(1,2,2)
%         plot(time_fr,nanmean(incNormFRs_var),'r')
%         hold on, grid on
%         plot(time_fr,nanmean(decNormFRs_var),'b')
%         plot(time_fr,nanmean(nonNormFRs_var),'k')
%         % plot(time_fr,mean(inhNormFRs_var),'g')
%         title('Variance in Normalized Firing Rates')
%         xlim([-2,0])
%         legend('Inc','Dec','Non','location','northwest')
%         saveas(gcf,[fpath,fname,'/NormalizedFRs'],'epsc')
% 
%         
%         try
%             figure
%             subplot(1,2,1)
%             plot_signal_ci(time_fr,(incNormFRs),'r')
%             hold on
%             plot_signal_ci(time_fr,(decNormFRs),'b')
%             hold on
%             plot_signal_ci(time_fr,(nonNormFRs),'k')
%             title('Normalized FRs + 95% CI')
%             xlim([-2,0])
%             legend('Inc','Dec','Non','location','northwest')
% 
% 
%             subplot(1,2,2)
%             plot_signal_ci(time_fr,(incNormFRs_var),'r')
%             hold on
%             plot_signal_ci(time_fr,(decNormFRs_var),'b')
%             hold on
%             plot_signal_ci(time_fr,(nonNormFRs_var),'k')
%             title('Normalized FR Variance + 95% CI')
%             xlim([-2,0])
%             legend('Inc','Dec','Non','location','northwest')
%             saveas(gcf,[fpath,fname,'/NormalizedFRs95ci'],'epsc')
%         catch
%             warning('Only one inc, dec, or non; cannot do CI');
%         end
%         
%         figure
%         
%         plot(time_fr,nanmean(incFRs),'r')
%         hold on, grid on
%         plot(time_fr,nanmean(decFRs),'b')
%         plot(time_fr,nanmean(nonFRs),'k')
%         title('Raw Firing Rates')
%         xlim([-2,0])
%         legend('Increasing','Decreasing','Non-Changing','location','northwest')
% 
%         saveas(gcf,[fpath,fname,'/FRs'],'epsc')
% 
% 
% 
%         % subplot(2,2,3)
%         % % map = brewermap(3,'Set1'); 
%         % histogram(incTaus(incTaus < 500),'facecolor','r','facealpha',.5,'edgecolor','none','Normalization','probability')
%         % hold on
%         % histogram(decTaus(decTaus < 500),'facecolor','b','facealpha',.5,'edgecolor','none','Normalization','probability')
%         % histogram(nonTaus(nonTaus < 500),'facecolor','k','facealpha',.5,'edgecolor','none','Normalization','probability')
%         % box off
%         % axis tight
%         % legend('Inc','Dec','Non','location','northeast')
% 
% 
%         % subplot(2,2,3)
% 
%         % max(length(decTaus(decTaus < 500)),length(incTaus(incTaus < 500))
%         tau_cutoff = 500;
%         nInc2 = length(incTaus(incTaus < tau_cutoff));
%         nDec2 = length(decTaus(decTaus < tau_cutoff));
%         nNon2 = length(nonTaus(nonTaus < tau_cutoff));
%         taumat = nan(max([nInc2,nDec2,nNon2]),3);
% 
% 
% 
%         taumat(1:nInc2,1) = incTaus(incTaus < tau_cutoff);
%         taumat(1:nDec2,2) = decTaus(decTaus < tau_cutoff);
%         taumat(1:nNon2,3) = nonTaus(nonTaus < tau_cutoff);
% 
% %         figure
% %         subplot(2,2,1)
% %         violin(taumat,'xlabel',{'inc','dec','nochange'},...
% %             'facecolor',[1 0 0;0 0 1;0 0 0]);
% %         % ylim([0 300])
% %         subplot(2,2,2)
% %         histogram(incTaus(incTaus < tau_cutoff)),title('Increasing Ns Taus')
% % 
% %         subplot(2,2,3)
% %         histogram(decTaus(decTaus<tau_cutoff)),title('Decreasing Ns Taus')
% % 
% %         subplot(2,2,4)
% %         histogram(nonTaus(nonTaus<tau_cutoff)),title('NoChange Ns Taus')
% %         saveas(gcf,[fpath,fname,'/TauFrequencies'],'epsc')
% 
%         try
%             figure
%             % Increasing
% 
%             tau_cutoff = 500;
%             R2_cutoff = 0.4;
%             p_cutoff = 0.01;
%             idx = (incTaus < tau_cutoff) & (incR2 > R2_cutoff) & (incP < p_cutoff);
%             % tbl = table(incTaus(idx),incGain(idx));
%             mdl = fitlm(incTaus(idx),incGain(idx));
%             subplot(1,2,1)
%             plot(mdl)
%             legend
%             xlabel('Autocorrelation Timescale'),ylabel('Fitted Gain Parameter')
%             title(['Tau vs Ramping for Increasing Neurons, ',num2str(sum(idx)),' valid Ns'])
%             % Decreasing
% 
%             idx = (decTaus < tau_cutoff) & (decR2 > R2_cutoff) & (decP < p_cutoff);
%             tbl = table(decTaus(idx),decGain(idx));
%             mdl = fitlm(decTaus(idx),decGain(idx));
%             subplot(1,2,2)
%             plot(mdl)
%             legend
%             xlabel('Autocorrelation Timescale'),ylabel('Fitted Gain Parameter')
%             title(['Tau vs Ramping for Decreasing Neurons, ',num2str(sum(idx)),' valid Ns'])
%             saveas(gcf,[fpath,fname,'/tau_vs_gain_500'],'epsc')
%         catch
%             warning('Problem w/ plotting')
%         end
% 
% 
%         %% Gain Split
%         figure('Position', [10 10 900 600])
% 
%         subplot(2,3,1)
%         histogram(incGain(incGain < 2*iqr(incGain)),'FaceColor','r')
%         title(['Gain for ',num2str(length(incGain)),' inc neurons']);
% 
%         subplot(2,3,2)
%         gain_quartiles = quantile(incGain,[0,0.25,0.5,0.75,1]);
%         for pp = 1:length(gain_quartiles)-1
%             gain_ix = (incGain > gain_quartiles(pp)) & (incGain < gain_quartiles(pp+1));
%             plot(time_fr,nanmean(incNormFRs(gain_ix,:)),'LineWidth',2)
%             hold on
%         end
%         xlabel('Time (s)')
%         xlim([time_fr(1) time_fr(end)])
%         ylabel('Normalized FR')
%         legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
%         title('Normalized FRs split via gain')
% 
%         subplot(2,3,3)
%         for pp = 1:length(gain_quartiles)-1
%             gain_ix = (incGain > gain_quartiles(pp)) & (incGain < gain_quartiles(pp+1));
%             plot(time_fr,nanmean(incFRs(gain_ix,:)),'LineWidth',2)
%             hold on
%         end
%         xlabel('Time (s)')
%         xlim([time_fr(1) time_fr(end)])
%         ylabel('FR (Hz)')
%         legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
%         title('Raw FRs split via gain')
% 
%         % Decreasing Neurons
%         subplot(2,3,4)
%         histogram(decGain(decGain < 2*iqr(decGain)),'FaceColor','b')
%         title(['Gain for ',num2str(length(decGain)),' dec neurons']);
% 
%         subplot(2,3,5)
%         gain_quartiles = quantile(decGain,[0,0.25,0.5,0.75,1]);
%         for pp = 1:length(gain_quartiles)-1
%             gain_ix = (decGain > gain_quartiles(pp)) & (decGain < gain_quartiles(pp+1));
%             plot(time_fr,nanmean(decNormFRs(gain_ix,:)),'LineWidth',2)
%             hold on
%         end
%         xlabel('Time (s)')
%         xlim([time_fr(1) time_fr(end)])
%         ylabel('Normalized FR')
%         legend({'Q1','Q2','Q3','Q4'},'Location','southwest','FontSize',14)
%         title('Normalized FRs split via gain')
% 
%         subplot(2,3,6)
%         for pp = 1:length(gain_quartiles)-1
%             gain_ix = (decGain > gain_quartiles(pp)) & (decGain < gain_quartiles(pp+1));
%             plot(time_fr,nanmean(decFRs(gain_ix,:)),'LineWidth',2)
%             hold on
%         end
%         xlabel('Time (s)')
%         xlim([time_fr(1) time_fr(end)])
%         ylabel('FR (Hz)')
%         legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
%         title('Raw FRs split via gain')
%         saveas(gcf,[fpath,fname,'/GainSplit'],'epsc')

        %% Plotting Spike Rasters


        close all
        fprintf(t_fid,[t_fname,'\n']);
    
    end
end
        

