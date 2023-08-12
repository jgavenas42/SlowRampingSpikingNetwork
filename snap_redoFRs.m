
%% Initial
clear all
close all

addpath(genpath('./utils'))

%% Set-up

nettype = 'flux_eegfix';
Ksyn2use = 0.5;


fpath = ['./data/' nettype '/Ksyn_' num2str(Ksyn2use) '/'];

files = dir([fpath,'*.mat']);


% Loading for pre-allocation
load([fpath files(1).name],'spks','time','norms')


for ii=1:length(files)
    
    tic
    fname = files(ii).name;
    load([fpath,fname],'spks','time','time_fr','Ksyn','network')
    
    
    type = cell(size(spks,2),1);
    gain = cell(size(spks,2),1);
    taus = nan(size(spks,2),1);
    acfs = nan(size(spks,2),12);
    
    norms = nan(size(spks,2),1);
    
    mean_frs = nan(size(spks,2),1);
    
    Zs = nan(size(spks,2),1);

    
    if size(spks,1) < 10 % If less than 10 trials...
%         type(:) = {'bad'};
%         incSpks = nan(3500,1);
%         decSpks = nan(3500,1);
%         nonSpks = nan(3500,1);
%         save([fpath,fname],'type','incSpks','decSpks','nonSpks','-append')
        continue
    end

    spks_sm = smoothdata(spks,3,'gaussian',400)*1000;
    
    acfs_all = zeros(size(spks,2),1001);
    
    for jj = 1:size(spks,2)
        spks_this = squeeze(spks(:,jj,:));
        
        idxs = ((time >= -3000) & (time < -2000));
        baseFRs = mean(spks_this(:,idxs),2)*1000;
        
        idxs = ((time >= -400) & (time < 0));
        stopFRs = nanmean(spks_this(:,idxs),2)*1000;
        

        % Rank-sum test to classify
        [P,H,STATS] = ranksum(stopFRs,baseFRs,'alpha',0.01);

        
        if isnan(H)
            type{jj} = 'bad';
        elseif ~H
            type{jj} = 'nochange';
        elseif (STATS.zval > 0)
            type{jj} = 'inc';
        elseif (STATS.zval < 0)
            type{jj} = 'dec';
        end
        
        Zs(jj) = STATS.zval;
        
        [taus(jj), acfs(jj,:)] = calculateTimescale('spk_mat',spks_this(:,501:1500));
        
        
        % Normalization to categorize heterogeneity
        spks_this_sm = squeeze(spks_sm(:,jj,1:3000));
        spks_this_sm = spks_this_sm(~isnan(spks_this_sm(:,1)),:);
        spks_this_sm = spks_this_sm / max(spks_this_sm,[],'all');
        
        spks_this_mean = squeeze(nanmean(spks_sm(:,jj,1:3000),1));
        
        nCk = nchoosek(1:size(spks_sm,1),2);

        norms_this = 0;
        for kk=1:length(nCk)
    %         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
            norms_this = norms_this + norm(spks_this_sm(nCk(kk,1))-spks_this_sm(nCk(kk,2)));

        end
%         norms_this = 0;
%         for kk=1:size(spks_this_sm,1)
%             norms_this = norms_this + norm(spks_this_sm(kk)-spks_this_mean);
%         end
        norms(jj) = norms_this / length(nCk);
        
        mean_frs(jj) = mean(spks_this_sm,'all');
        
        spks_this_sm = smoothdata(spks_this,2,'movmean',10);
        
        spks_this_sm = spks_this_sm(~isnan(spks_this_sm(:,1)),:);
        
        acfs_this = zeros(size(spks_this_sm,1),1001);
        for kk = 1:size(spks_this_sm,1)
            acfs_this(kk,:) = autocorr(spks_this_sm(kk,:)',1000);
        end
        acfs_all(jj,:) = nanmean(acfs_this,1);
    end
    
    
    spks_sm = smoothdata(spks,3,'movmean',10);
    cluster_FRs = squeeze(mean(reshape(spks_sm,size(spks,1),80,4,[]),2));
    
    cluster_FRs = cluster_FRs(~isnan(cluster_FRs(:,1,1)),:,:);
    acfs_clust = zeros(4,1001);
    for jj = 1:size(cluster_FRs,2)
        frs_this = squeeze(cluster_FRs(:,jj,:));
        acfs_this = zeros(size(frs_this,1),1001);
        for kk = 1:size(frs_this,1)
            acfs_this(kk,:) = autocorr(frs_this(kk,:)',1000);
        end
        acfs_clust(jj,:) = nanmean(acfs_this,1);
    end
    
    spks_sm = smoothdata(spks,3,'gaussian',400)*1000;
    nInc = sum(strcmp(type,'inc'));
    nDec = sum(strcmp(type,'dec'));
    nNon = sum(strcmp(type,'non'));
    
    incSpks = squeeze(nanmean(spks_sm(:,strcmp(type,'inc'),:),[1 2]));
    decSpks = squeeze(nanmean(spks_sm(:,strcmp(type,'dec'),:),[1 2]));
    nonSpks = squeeze(nanmean(spks_sm(:,strcmp(type,'non'),:),[1 2]));
    
    acf = nanmean(acfs,1);
    
    
    toc
    
%     save([fpath,fname],'type','incSpks','decSpks',...
%         'nonSpks','taus','acfs','norms','mean_frs',...
%         'acfs_all','acfs_clust','-append')
%     save([fpath,fname],'type','incSpks','decSpks','nonSpks','-append')
    
    
    [t_, inc_ix] = min(abs(Zs - nanmean(Zs(Zs > 1.96))));
    [t_, dec_ix] = min(abs(Zs - nanmean(Zs(Zs < -1.96))));
    
    close all
    
%     figure('Position',[0 0 300 300])
%     
%     tiledlayout(2,1,'TileSpacing','compact')
%     nexttile
%     plotRaster(squeeze(spks(~isnan(spks(:,1,1)),inc_ix,:)),'Color',[1 0.3 0.2],'MarkerSize',5)
%     ylim([0,size(spks(~isnan(spks(:,1,1)),:,:),1)+1])
%     xlim([0 3000])
%     ax = gca;
%     axis off
%     ax.YAxis.Visible = 'on';
%     ylabel('Trial','FontSize',16)
%     
%     nexttile
%     plot(time/1000,squeeze(nanmean(spks_sm(:,inc_ix,:),1)),'Color',[1 0.3 0.2],'LineWidth',3)
%     hold on
% 
%     xlim([-3 0.])
%     ax = gca;
%     axis off
%     ax.YAxis.Visible = 'on';
%     ax.XAxis.Visible = 'on';
%     ylabel('Firing Rate (Hz)','FontSize',16)
%     xlabel('Time rel. threshold-crossing (s)','FontSize',16)
%     
%     saveas(gcf,['./figures/examples/' nettype '_RasterInc_Ksyn_',num2str(Ksyn), '_Net_',num2str(network),'.png'])


%     figure('Position',[0 0 300 300])
%     
%     tiledlayout(2,1,'TileSpacing','compact')
%     nexttile
%     plotRaster(squeeze(spks(~isnan(spks(:,1,1)),dec_ix,:)),'Color',[0 0.5 1],'MarkerSize',5)
%     ylim([0,size(spks(~isnan(spks(:,1,1)),:,:),1)+1])
%     xlim([0 3000])
%     ax = gca;
%     axis off
%     ax.YAxis.Visible = 'on';
%     ylabel('Trial','FontSize',16)
%     
%     nexttile
%     plot(time/1000,squeeze(nanmean(spks_sm(:,dec_ix,:),1)),'Color',[0 0.5 1],'LineWidth',3)
%     hold on
%     
%     xlim([-3 0.])
%     ax = gca;
%     axis off
%     ax.YAxis.Visible = 'on';
%     ax.XAxis.Visible = 'on';
%     ylabel('Firing Rate (Hz)','FontSize',16)
%     xlabel('Time rel. threshold-crossing (s)','FontSize',16)
%     
%     
%     
%     saveas(gcf,['./figures/examples/' nettype '_RasterDec_Ksyn_',num2str(Ksyn), '_Net_',num2str(network),'.png'])
    
    
    
    
end