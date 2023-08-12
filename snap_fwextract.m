tic
clear all
addpath(genpath('utils/'))
fpath = ('data_real/');

filenames = dir([fpath,'*.mat']);

load([fpath,filenames(1).name]);

figures_path = 'figures/scratch';

% Return all of the firing rates in ___FRs_full variable?
return_full = false;

baseline=500;

nInc = 0;
nDec = 0;
nNon = 0;

for i=1:length(filenames)
    load([fpath,filenames(i).name],'type')
    
    if strcmp(type,'inc')
        nInc = nInc + 1;
    elseif strcmp(type,'dec')
        nDec = nDec + 1;
    elseif strcmp(type,'nochange')
        nNon = nNon + 1;
    end
end



load([fpath,filenames(1).name]);

time_fr = linspace(-3.5,0,size(normFRs,2));


decFRs = zeros(nDec,size(normFRs,2));
incFRs = zeros(nInc,size(normFRs,2));
nonFRs = zeros(nNon,size(normFRs,2));

incSpks = zeros(nInc,3501);
decSpks = zeros(nDec,3501);
time = linspace(-3.5,0,3501);

decNormFRs = zeros(nDec,size(normFRs,2));
incNormFRs = zeros(nInc,size(normFRs,2));
nonNormFRs = zeros(nNon,size(normFRs,2));

decNormFRs_var = zeros(nDec,size(normFRs,2));
incNormFRs_var = zeros(nInc,size(normFRs,2));
nonNormFRs_var = zeros(nNon,size(normFRs,2));

decTaus = zeros(nDec,1);
incTaus = zeros(nInc,1);
nonTaus = zeros(nNon,1);

decRegion = cell(nDec,1);
incRegion = cell(nInc,1);
nonRegion = cell(nNon,1);

decIter = 0;
incIter = 0;
nonIter = 0;

incGain = zeros(nInc,1);
decGain = zeros(nDec,1);

incR2 = zeros(nInc,1);
decR2 = zeros(nDec,1);

incP = zeros(nInc,1);
decP = zeros(nDec,1);

decBI = zeros(nDec,1);
incBI = zeros(nInc,1);
nonBI = zeros(nNon,1);

norms_inc_real = nan(nInc,1);
norms_dec_real = nan(nDec,1);

decACF = zeros(nDec,1001);
incACF = zeros(nInc,1001);

time_decim = time(500:50:end);

incSpks_sig = zeros(nInc,length(time_decim),2);
decSpks_sig = zeros(nDec,length(time_decim),2);

% incFitted = cell{nInc,2}; % First column fit objects, second gof object
% decFitted = cell{nDec,2}; % First column fit objects, second column gof object

if return_full
    incFRs_full = [];
    decFRs_full = [];
end

for i = 1:length(filenames)
    
    load([fpath,filenames(i).name],'normFRs','type','BaselineAlignedSpikes','p','tau','FRs','region','BI','avgdist','StopAlignedSpikes','baseFRs')
    
    if strcmp(type,'bad'),continue,end
    
    tmpFR = nanmean(normFRs);
    
    tmpFRvar = nanvar(normFRs);
    
    tmpFRfull = nanmean(FRs);
    
%     tmpTau = calculateTimescale('spk_mat',BaselineAlignedSpikes);
    tmpTau = tau;
    tmpBI = BI;
    
    
    
    if strcmp(type,'inc')
        incIter = incIter + 1;
        
        incFRs(incIter,:) = tmpFRfull;
        
        spks_sm = smoothdata(StopAlignedSpikes,2,'gaussian',400) * 400;
        incSpks(incIter,:) = nanmean(spks_sm,1) - nanmean(spks_sm(:,1:baseline),'all');
        
        incNormFRs(incIter,:) = tmpFR;
        incNormFRs_var(incIter,:) = tmpFRvar;
        incTaus(incIter) = tmpTau;
        
        incBI(incIter) = tmpBI;
        
        incP(incIter) = p;
        
        incRegion{incIter} = region;
        
        load([fpath,filenames(i).name],'gain','sigmoid_gof','avgdist')
        incGain(incIter) = gain;
        incR2(incIter) = sigmoid_gof.rsquare;
        
        norms_inc_real(incIter) = avgdist;
        
        if return_full
%             if gain > 0.2858
%             if strcmp(region,'LpSMA')
                incFRs_full = [incFRs_full;FRs];
%             end
        end
        
        spks_this_sm = smoothdata(StopAlignedSpikes,2,'movmean',10);
        
        spks_this_sm = spks_this_sm(~isnan(spks_this_sm(:,1)),:);
        
        acfs_this = zeros(size(spks_this_sm,1),1001);
        for kk = 1:size(spks_this_sm,1)
            acfs_this(kk,:) = autocorr(spks_this_sm(kk,:)',1000);
        end
        incACF(incIter,:) = nanmean(acfs_this,1);
        
        
        
        for kk = 1:length(time_decim)
            ix2use = find(time > time_decim(kk)-1 & time <= time_decim(kk));
            lm = fitlm(time(ix2use),incSpks(incIter,ix2use));
            incSpks_sig(incIter,kk,1) = lm.Coefficients.pValue(2);
            incSpks_sig(incIter,kk,2) = lm.Coefficients.Estimate(2);
        end
        
%         close all
%         figure('Position',[0 0 300 300])
%         tiledlayout(2,1,'TileSpacing','compact')
%         nexttile
%         plotRaster(squeeze(StopAlignedSpikes),'Color',[1 0.3 0.2],'MarkerSize',5)
%         ylim([0,size(StopAlignedSpikes,1)+1])
%         xlim([501 3500])
%         ax = gca;
%         axis off
%         ax.YAxis.Visible = 'on';
%         ylabel('Trial','FontSize',16)
% 
%         nexttile
%         plot(time,squeeze(nanmean(spks_sm,1)),'Color',[1 0.3 0.2],'LineWidth',3)
%         hold on
% 
%         xlim([-3 0])
%         ax = gca;
%         axis off
%         ax.YAxis.Visible = 'on';
%         ax.XAxis.Visible = 'on';
%         ylabel('Firing Rate (Hz)','FontSize',16)
%         xlabel('Time rel. movement (s)','FontSize',16)
%         saveas(gcf,['./figures/examples_real/Inc_Raster_',num2str(incIter), '.png'])
        
    elseif strcmp(type,'dec')
        decIter = decIter + 1;
        
        decFRs(decIter,:) = tmpFRfull;
        
        spks_sm = smoothdata(StopAlignedSpikes,2,'gaussian',400) * 400;
        decSpks(decIter,:) = nanmean(spks_sm,1) - nanmean(spks_sm(:,1:baseline),'all');
        
        decNormFRs(decIter,:) = tmpFR;
        decNormFRs_var(decIter,:) = tmpFRvar;
        decTaus(decIter,:) = tmpTau;
        
        decBI(decIter) = tmpBI;
        
        decP(decIter) = p;
        
        decRegion{decIter} = region;
        
        load([fpath,filenames(i).name],'gain','sigmoid_gof','avgdist')
        decGain(decIter) = gain;
        decR2(decIter) = sigmoid_gof.rsquare;
        
        norms_dec_real(decIter) = avgdist;
        
        if return_full
%             if gain > 0.1994
%             if strcmp(region,'LpSMA')
                decFRs_full = [decFRs_full;FRs];
%             end
        end
        
        spks_this_sm = smoothdata(StopAlignedSpikes,2,'movmean',10);
        
        spks_this_sm = spks_this_sm(~isnan(spks_this_sm(:,1)),:);
        
        acfs_this = zeros(size(spks_this_sm,1),1001);
        for kk = 1:size(spks_this_sm,1)
            acfs_this(kk,:) = autocorr(spks_this_sm(kk,:)',1000);
        end
        decACF(decIter,:) = nanmean(acfs_this,1);
        
        for kk = 1:length(time_decim)
            ix2use = find(time > time_decim(kk)-1 & time <= time_decim(kk));
            lm = fitlm(time(ix2use),decSpks(decIter,ix2use));
            decSpks_sig(decIter,kk,1) = lm.Coefficients.pValue(2);
            decSpks_sig(decIter,kk,2) = lm.Coefficients.Estimate(2);
        end
        
%         close all
%         figure('Position',[0 0 300 300])
%         tiledlayout(2,1,'TileSpacing','compact')
%         nexttile
%         plotRaster(squeeze(StopAlignedSpikes),'Color',[0 0.5 1],'MarkerSize',5)
%         ylim([0,size(StopAlignedSpikes,1)+1])
%         xlim([501 3500])
%         ax = gca;
%         axis off
%         ax.YAxis.Visible = 'on';
%         ylabel('Trial','FontSize',16) 
% 
%         nexttile
%         plot(time,squeeze(nanmean(spks_sm,1)),'Color',[0 0.5 1],'LineWidth',3)
%         hold on
% 
%         xlim([-3. 0.])
%         ax = gca;
%         axis off
%         ax.YAxis.Visible = 'on';
%         ax.XAxis.Visible = 'on';
%         ylabel('Firing Rate (Hz)','FontSize',16)
%         xlabel('Time rel. movement (s)','FontSize',16)
%         saveas(gcf,['./figures/examples_real/Dec_Raster_',num2str(decIter), '.png'])

        
    elseif strcmp(type,'nochange')
        nonIter = nonIter + 1;
        
        nonFRs(nonIter,:) = tmpFRfull;
        
        nonBI(nonIter) = tmpBI;
        
        nonRegion{nonIter} = region;
        
        nonNormFRs(nonIter,:) = tmpFR;
        nonNormFRs_var(nonIter,:) = tmpFRvar;
        nonTaus(nonIter,:) = tmpTau;
        
    end
end

if return_full
    incFRs_full = smoothdata(incFRs_full,2,'gaussian',8);
    decFRs_full = smoothdata(decFRs_full,2,'gaussian',8);
end
toc

mfcRegions = unique([incRegion; decRegion; nonRegion]);

mfcRegions = mfcRegions([2 3 4 6 12 13 16 17 18 19 28 29 30]);

disp(['# MFC Neurons: ' num2str(sum(ismember(incRegion,mfcRegions)) + sum(ismember(decRegion,mfcRegions)) + sum(ismember(nonRegion,mfcRegions)))]);
disp([num2str(sum(ismember(incRegion,mfcRegions))), ' increasing, ', num2str(sum(ismember(decRegion,mfcRegions))), ' decreasing']);


% save('./heterogeneity/l2norms.mat',...
%     'norms_inc_real','norms_dec_real','-append')

%% Plotting Firing Rates
time_fr = -3.45:0.05:0;
figure

subplot(1,2,1)
plot(time_fr,nanmean(incNormFRs),'r')
hold on, grid on
plot(time_fr,nanmean(decNormFRs),'b')
plot(time_fr,nanmean(nonNormFRs),'k')
title('Normalized Firing Rates')
xlim([-2,0])
legend('Increasing','Decreasing','Non-Changing','location','northwest')


subplot(1,2,2)
plot(time_fr,nanmean(incNormFRs_var),'r')
hold on, grid on
plot(time_fr,nanmean(decNormFRs_var),'b')
plot(time_fr,nanmean(nonNormFRs_var),'k')
title('Variance in Normalized Firing Rates')
xlim([-2,0])
legend('Increasing','Decreasing','Non-Changing','location','northwest')
saveas(gcf,[figures_path,'/normalizedFRs'],'epsc')

figure
subplot(1,2,1)
plot_signal_ci(time_fr,(incNormFRs),'r')
hold on
plot_signal_ci(time_fr,(decNormFRs),'b')
hold on
plot_signal_ci(time_fr,(nonNormFRs),'k')
title('Normalized FRs + 95% CI')
xlim([-2,0])
legend('Increasing','Decreasing','Non-Changing','location','northwest')


subplot(1,2,2)
plot_signal_ci(time_fr,(incNormFRs_var),'r')
hold on
plot_signal_ci(time_fr,(decNormFRs_var),'b')
hold on
plot_signal_ci(time_fr,(nonNormFRs_var),'k')
title('Normalized FR Variance + 95% CI')
xlim([-2,0])
legend('Increasing','Decreasing','Non-Changing','location','northwest')
saveas(gcf,[figures_path,'/normalizedFRs95ci'],'epsc')


figure
plot(time_fr,nanmean(incFRs),'r')
hold on, grid on
plot(time_fr,nanmean(decFRs),'b')
plot(time_fr,nanmean(nonFRs),'k')
title('Raw Firing Rates')
xlim([-2,0])
legend('Increasing','Decreasing','Non-Changing','location','northwest')

% saveas(gcf,[figures_path,'/FRs'],'epsc')



% subplot(2,2,3)
% % map = brewermap(3,'Set1'); 
% histogram(incTaus(incTaus < 500),'facecolor','r','facealpha',.5,'edgecolor','none','Normalization','probability')
% hold on
% histogram(decTaus(decTaus < 500),'facecolor','b','facealpha',.5,'edgecolor','none','Normalization','probability')
% histogram(nonTaus(nonTaus < 500),'facecolor','k','facealpha',.5,'edgecolor','none','Normalization','probability')
% box off
% axis tight
% legend('Inc','Dec','Non','location','northeast')


% subplot(2,2,3)

figure
% max(length(decTaus(decTaus < 500)),length(incTaus(incTaus < 500))
tau_cutoff = 500;
nInc2 = length(incTaus(incTaus < tau_cutoff));
nDec2 = length(decTaus(decTaus < tau_cutoff));
nNon2 = length(nonTaus(nonTaus < tau_cutoff));
taumat = nan(max([nInc2,nDec2,nNon2]),3);



taumat(1:nInc2,1) = incTaus(incTaus < tau_cutoff);
taumat(1:nDec2,2) = decTaus(decTaus < tau_cutoff);
taumat(1:nNon2,3) = nonTaus(nonTaus < tau_cutoff);

subplot(2,2,1)
violin(taumat,'xlabel',{'inc','dec','nochange'},...
    'facecolor',[1 0 0;0 0 1;0 0 0]);
% ylim([0 300])
subplot(2,2,2)
histogram(incTaus(incTaus < tau_cutoff)),title('Increasing Ns Taus')

subplot(2,2,3)
histogram(decTaus(decTaus<tau_cutoff)),title('Decreasing Ns Taus')

subplot(2,2,4)
histogram(nonTaus(nonTaus<tau_cutoff)),title('NoChange Ns Taus')

% saveas(gcf,[figures_path,'/TauFrequencies'],'epsc')

%% Splitting firing rates based on gain

figure('Position', [10 10 900 600])

subplot(2,3,1)
histogram(incGain(incGain < 2*iqr(incGain)),'FaceColor','r')
title(['Gain for ',num2str(length(incGain)),' inc neurons']);

subplot(2,3,2)
gain_quartiles = quantile(incGain,[0,0.25,0.5,0.75,1]);
for pp = 1:length(gain_quartiles)-1
    gain_ix = (incGain > gain_quartiles(pp)) & (incGain < gain_quartiles(pp+1));
    plot(time_fr,nanmean(incNormFRs(gain_ix,:)),'LineWidth',2)
    hold on
end
xlabel('Time (s)')
xlim([time_fr(1) time_fr(end)])
ylabel('Normalized FR')
legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
title('Normalized FRs split via gain')

% subplot(2,3,3)
% for pp = 1:length(gain_quartiles)-1
%     gain_ix = (incGain > gain_quartiles(pp)) & (incGain < gain_quartiles(pp+1));
%     plot(time_fr,nanmean(incFRs(gain_ix,:)),'LineWidth',2)
%     hold on
% end
% xlabel('Time (s)')
% xlim([time_fr(1) time_fr(end)])
% ylabel('FR (Hz)')
% legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
% title('Raw FRs split via gain')

subplot(2,3,3)
for pp = 1:length(gain_quartiles)-1
    gain_ix = (incGain > gain_quartiles(pp)) & (incGain < gain_quartiles(pp+1));
    plot(time_fr,nanmean(incNormFRs_var(gain_ix,:)),'LineWidth',2)
    hold on
end
xlabel('Time (s)')
xlim([time_fr(1) time_fr(end)])
ylabel('Variance')
legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
title('Variance of Normalized Firing Rates')

% Decreasing Neurons
subplot(2,3,4)
histogram(decGain(decGain < 2*iqr(decGain)),'FaceColor','b')
title(['Gain for ',num2str(length(decGain)),' dec neurons']);

subplot(2,3,5)
gain_quartiles = quantile(decGain,[0,0.25,0.5,0.75,1]);
for pp = 1:length(gain_quartiles)-1
    gain_ix = (decGain > gain_quartiles(pp)) & (decGain < gain_quartiles(pp+1));
    plot(time_fr,nanmean(decNormFRs(gain_ix,:)),'LineWidth',2)
    hold on
end
xlabel('Time (s)')
xlim([time_fr(1) time_fr(end)])
ylabel('Normalized FR')
legend({'Q1','Q2','Q3','Q4'},'Location','southwest','FontSize',14)
title('Normalized FRs split via gain')

% subplot(2,3,6)
% for pp = 1:length(gain_quartiles)-1
%     gain_ix = (decGain > gain_quartiles(pp)) & (decGain < gain_quartiles(pp+1));
%     plot(time_fr,nanmean(decFRs(gain_ix,:)),'LineWidth',2)
%     hold on
% end
% xlabel('Time (s)')
% xlim([time_fr(1) time_fr(end)])
% ylabel('FR (Hz)')
% legend({'Q1','Q2','Q3','Q4'},'Location','southwest','FontSize',14)
% title('Raw FRs split via gain')
subplot(2,3,6)
for pp = 1:length(gain_quartiles)-1
    gain_ix = (decGain > gain_quartiles(pp)) & (decGain < gain_quartiles(pp+1));
    plot(time_fr,nanmean(decNormFRs_var(gain_ix,:)),'LineWidth',2)
    hold on
end
xlabel('Time (s)')
xlim([time_fr(1) time_fr(end)])
ylabel('Variance')
legend({'Q1','Q2','Q3','Q4'},'Location','northwest','FontSize',14)
title('Variance of Normalized Firing Rates')

% saveas(gcf,[figures_path,'/GainSplit'],'epsc')



%% Plotting individual activities

figure
subplot(2,2,1)
plotSpikeRaster(WAlignedSpikes==1,'PlotType','vertline');

subplot(2,2,2)
time_fr = linspace(-3.5,0,size(normFRs,2));
plot(time_fr,firingRate(WAlignedSpikes)),hold on
plot(time_fr,nanmean(firingRate(WAlignedSpikes),1),'k','linewidth',2)
xlim([-2.5 0])

subplot(2,2,3)
violin([baseFRs,stopFRs],'xlabel',{'baseline FR','near-W FR'});
title('Firing Rate violin')

subplot(2,2,4)
psth(find(StopAlignedSpikes==1)',100,1000,ntrials,3501)


%% Scratch work for fitting sigmoids

% Individual trials
trial_num = 6;
fB = num2str(mean(WAlignedSpikes(trial_num,find((time > -2.5) & (time <= -1.5)))));
fW = num2str(mean(WAlignedSpikes(trial_num,find((time > -0.4)))));

model = strcat('(',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);

idxs = find((time > -2.5));
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Algorithm','Levenberg-Marquardt');
ft = fittype(model,'options',fo);

[fitted,gof,output] = fit(time(idxs)',smoothdata(WAlignedSpikes(trial_num,idxs),'gaussian',200)',ft);

plot(fitted,time(idxs),smoothdata(WAlignedSpikes(trial_num,idxs),'gaussian',200))
xlim([-2.5 0])
legend


%% All trials - fitting on firing rate

% good increasing idx = 
fB = num2str(mean(baseFRs));
fW = num2str(mean(stopFRs));

if strcmp(type,'inc')
    model = strcat('(',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
else
    model = strcat('(-',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
end

idxs = find((time_fr > -2.5));
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Algorithm','Levenberg-Marquardt');
ft = fittype(model,'options',fo);

frs = nanmean(firingRate(WAlignedSpikes),1);
[fitted,gof,output] = fit(time_fr(idxs)',smoothdata(frs(idxs),'gaussian',4)',ft);
figure
plot(fitted,time_fr(idxs),smoothdata(frs(idxs),'gaussian',4))
xlim([-2.5 0])
legend

%% Correlating the fitted gain (speed of rise) against autocorrelation timescale
figure
% Increasing

tau_cutoff = 500;
R2_cutoff = 0.4;
p_cutoff = 0.01;
idx = (incTaus < tau_cutoff) & (incR2 > R2_cutoff) & (incP < p_cutoff);
% tbl = table(incTaus(idx),incGain(idx));
mdl = fitlm(incTaus(idx),incGain(idx));
subplot(1,2,1)
plot(mdl)
legend
xlabel('Autocorrelation Timescale'),ylabel('Fitted Gain Parameter')
title(['Tau vs Ramping for Increasing Neurons, ',num2str(sum(idx)),' valid Ns'])
% Decreasing

idx = (decTaus < tau_cutoff) & (decR2 > R2_cutoff) & (decP < p_cutoff);
tbl = table(decTaus(idx),decGain(idx));
mdl = fitlm(decTaus(idx),decGain(idx));
subplot(1,2,2)
plot(mdl)
legend
xlabel('Autocorrelation Timescale'),ylabel('Fitted Gain Parameter')
title(['Tau vs Ramping for Decreasing Neurons, ',num2str(sum(idx)),' valid Ns'])

saveas(gcf,[figures_path,'/tauVSgain_',num2str(tau_cutoff)],'epsc')



%% Compare regions for increasing neurons

smooth_window = 4;
[unIncRegions,ia,ic] = unique(incRegion);
unIncRegionCounts = accumarray(ic,1);

figure('Position', [10 10 900 400])

sgtitle('Inc Neurons: Frequency & Single-Neuron Activity')

subplot(2,3,1)
bar(categorical(unIncRegions),unIncRegionCounts)
title('# Increasing Neurons From Brain Regions')
xlabel('Region');ylabel('# Neurons')

chosen_regions = unIncRegions(unIncRegionCounts >= 7);

subplot(2,3,2)
plot(time_fr,smoothdata(incNormFRs(strcmp(incRegion,'LSMA'),:),2,'movmean',smooth_window))
title('LSMA (N=10)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,3,3)
plot(time_fr,smoothdata(incNormFRs(strcmp(incRegion,'LpSMA'),:),2,'movmean',smooth_window))
title('LpSMA (N=19)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,3,4)
plot(time_fr,smoothdata(incNormFRs(strcmp(incRegion,'RSMA'),:),2,'movmean',smooth_window))
title('RSMA (N=11)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,3,5)
plot(time_fr,smoothdata(incNormFRs(strcmp(incRegion,'LAC'),:),2,'movmean',smooth_window))
title('LAC (N=7)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,3,6)
plot(time_fr,smoothdata(incNormFRs(strcmp(incRegion,'REC'),:),2,'movmean',smooth_window))
title('REC (N=7)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')
saveas(gcf,[figures_path,'/Regions_inc_activity'],'epsc')


%% Compare regions for decreasing neurons
[unDecRegions,ia,ic] = unique(decRegion);
unDecRegionCounts = accumarray(ic,1);

chosen_regions = unDecRegions(unDecRegionCounts >= 17);


figure('Position', [10 10 1200 400])

sgtitle('Dec Neurons: Frequency & Single-Neuron Activity')

subplot(2,4,1)
bar(categorical(unDecRegions),unDecRegionCounts)
title('# Decreasing Neurons From Brain Regions')
xlabel('Region');ylabel('# Neurons')

subplot(2,4,2)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'LpSMA'),:),2,'movmean',smooth_window))
title('LpSMA (N=50)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,3)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'RMC'),:),2,'movmean',smooth_window))
title('RMC (N=24)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,4)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'RSMA'),:),2,'movmean',smooth_window))
title('RSMA (N=26)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,5)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'RMC'),:),2,'movmean',smooth_window))
title('RMC (N=21)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,6)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'LOF'),:),2,'movmean',smooth_window))
title('LOF (N=21)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,7)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'ROF'),:),2,'movmean',smooth_window))
title('ROF (N=21)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')

subplot(2,4,8)
plot(time_fr,smoothdata(decNormFRs(strcmp(decRegion,'LAC'),:),2,'movmean',smooth_window))
title('LAC (N=17)')
xlim([-3.5 0])
xlabel('Time before movement (s)');ylabel('Firing Rate (Hz, smoothed)')
saveas(gcf,[figures_path,'/Regions_dec_activity'],'epsc')


%% Compare regions for non-changing neurons


figure('Position', [10 10 1200 400])

sgtitle('Distributions of changing vs non-changing neurons')

subplot(1,3,1)
[unNonRegions,ia,ic] = unique(nonRegion);
unNonRegionCounts = accumarray(ic,1);
bar(categorical(unNonRegions),unNonRegionCounts)
title(['Non-Changing Neurons (N=',num2str(nNon),')'])
xlabel('Region');ylabel('# Neurons')

subplot(1,3,2)
[unIncRegions,ia,ic] = unique(incRegion);
unIncRegionCounts = accumarray(ic,1);
bar(categorical(unIncRegions),unIncRegionCounts)
title(['Increasing Neurons (N=',num2str(nInc),')'])
xlabel('Region');ylabel('# Neurons')

subplot(1,3,3)
[unDecRegions,ia,ic] = unique(decRegion);
unDecRegionCounts = accumarray(ic,1);
bar(categorical(unDecRegions),unDecRegionCounts)
title(['Decreasing Neurons (N=',num2str(nDec),')'])
xlabel('Region');ylabel('# Neurons')

saveas(gcf,[figures_path,'/Regions'],'epsc')


%% Boxplot for Burst Indexes

figure('Position',[10 10 600 200])
% max(length(decTaus(decTaus < 500)),length(incTaus(incTaus < 500))

BImat = nan(max([nInc,nDec,nNon]),3);

ylims = [0 0.3];


BImat(1:nInc,1) = incBI;
BImat(1:nDec,2) = decBI;
BImat(1:nNon,3) = nonBI;

% subplot(1,3,1)
% violin(BImat,'xlabel',{'inc','dec','nochange'},...
%     'facecolor',[1 0 0;0 0 1;0 0 0]);
% ylim([0 300])
subplot(1,3,1)
histogram(decBI,'BinEdges',0:0.05:1,'FaceColor',[0 0 1],'FaceAlpha',0.5,'Normalization','probability')
title('Dec. Burst Index')
xline(nanmean(decBI));
text(nanmean(decBI)+0.01,0.25,['Mean = ',num2str(nanmean(decBI))]);
ylim(ylims)

subplot(1,3,2)
histogram(nonBI,'BinEdges',0:0.05:1,'FaceColor',[0 0 0],'FaceAlpha',0.5,'Normalization','probability')
title('NoChange Ns Burst Index')
xline(nanmean(nonBI));
text(nanmean(nonBI)+0.01,0.25,['Mean = ',num2str(nanmean(nonBI))]);
ylim(ylims)


subplot(1,3,3)
histogram(incBI,'BinEdges',0:0.05:1,'FaceColor',[1 0 0],'FaceAlpha',0.5,'Normalization','probability')
title('Inc. Burst Index')
xline(nanmean(incBI));
text(nanmean(incBI)+0.01,0.25,['Mean = ',num2str(nanmean(incBI))])
ylim(ylims)




% saveas(gcf,[figures_path,'/BIfreq'],'epsc')


%% Burst index vs gain parameter

% subplot(1,2,1)
% plot(decBI,decGain,'.b')
% title('Dec. Burst Index vs Gain')
% % xline(nanmean(decBI));
% % text(nanmean(decBI)+0.01,0.25,['Mean = ',num2str(nanmean(decBI))]);
% 
% % subplot(1,3,2)
% % plot(nonBI,nonGain,'.k')
% % title('No Change Burst Index vs Gain')
% % % xline(nanmean(nonBI));
% % % text(nanmean(nonBI)+0.01,0.25,['Mean = ',num2str(nanmean(nonBI))]);
% % % ylim(ylims)
% 
% 
% subplot(1,2,2)
% histogram(incBI,incGain,'.r')
% title('Inc. Burst Index vs Gain')
% xline(nanmean(incBI));
% text(nanmean(incBI)+0.01,0.25,['Mean = ',num2str(nanmean(incBI))])
% ylim(ylims)
figure('Position',[10 10 600 200])

gain_cutoff_high = 1.5;
gain_cutoff_low = 0.;
p_cutoff = 0.01;
idx = (incGain < gain_cutoff_high) & (incGain > gain_cutoff_low) & (incP < p_cutoff);
% tbl = table(incTaus(idx),incGain(idx));
subplot(1,2,1)
plot(incBI(idx),incGain(idx),'.r')
mdl = fitlm(incBI(idx),incGain(idx));
refline(mdl.Coefficients.Estimate)
legend({'data','linear model'})
xlabel('Burst Index'),ylabel('Fitted Gain Parameter')
title(['BI vs Gain for Increasing Neurons, ',num2str(sum(idx)),' valid Ns'])
text(0.4,0.6,['Fit R^2 = ',num2str(mdl.Rsquared.ordinary)])

% Decreasing

idx = (decGain < gain_cutoff_high) & (decGain > gain_cutoff_low) & (decP < p_cutoff);
subplot(1,2,2)
plot(decBI(idx),decGain(idx),'.b')
mdl = fitlm(decBI(idx),decGain(idx));
refline(mdl.Coefficients.Estimate)
legend({'data','linear model'})
xlabel('Burst Index'),ylabel('Fitted Gain Parameter')
title(['BI vs Gain for Decreasing Neurons, ',num2str(sum(idx)),' valid Ns'])
text(0.4,1,['Fit R^2 = ',num2str(mdl.Rsquared.ordinary)])


%% KS Density Plots - average
figure('Position', [10 10 750 500])
% ix_inc = incGain > prctile(incGain,66);
% ix_dec = decGain > prctile(decGain,66);
ix_inc = incGain > 0;
ix_dec = decGain > 0;


colors = [linspace(0.929,0.494,50);linspace(0.694,0.184,50);linspace(0.1250,0.556,50)]';
[f_base_inc,xi_base_inc] = ksdensity(nanmean(incFRs(ix_inc,21:40),2),0:0.5:50,'Support',[-0.001,100],'BoundaryCorrection','reflection');
[f_base_dec,xi_base_dec] = ksdensity(nanmean(decFRs(ix_inc,21:40),2),0:0.5:50,'Support',[-0.001,100],'BoundaryCorrection','reflection');

kldivs_inc = zeros(50,1);
kldivs_dec = zeros(50,1);
kldivs_cf = zeros(50,1);

colormap(colors)
for i=21:70
    subplot(2,3,1)
    [finc,xi] = ksdensity(incFRs(:,i),0:0.5:50,'Support',[-0.001,100],'BoundaryCorrection','reflection');
    kldivs_inc(i-20) = kldivergence(finc,f_base_inc);

    plot(xi,finc,'color',colors(i-20,:))
    hold on
    
    subplot(2,3,2)
    [fdec,xi]=ksdensity(decFRs(:,i),0:0.5:50,'Support',[-0.001,100],'BoundaryCorrection','reflection');
    kldivs_dec(i-20) = kldivergence(fdec,f_base_dec);
    plot(xi,fdec,'color',colors(i-20,:))
    hold on
    
    subplot(2,3,3)
    plot(xi,finc-fdec,'color',colors(i-20,:)), hold on
    kldivs_cf(i-20) = 0.5 * kldivergence(finc,fdec) + 0.5 * kldivergence(fdec,finc);
end

sgtitle('Density & Divergence of trial-averaged FRs')
subplot(2,3,1)
title("Density of Inc Neuron FRs over time")
xlabel('Firing Rate'),xlim([0 40])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');
subplot(2,3,2)
title("Density of Dec Neuron FRs over time")
xlabel('Firing Rate'),xlim([0 20])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');

subplot(2,3,3)
title("Inc - Dec Neuron Densitities over time")
xlabel('Firing Rate'),xlim([0 20])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');

subplot(2,3,4)
title("Inc: single-time FRs vs baseline")
plot(time_fr(time_fr > -2.5),kldivs_inc)
yl = get(gca,'YLim');
patch([-2.5 -1.5 -1.5 -2.5],[yl(1) yl(1) yl(2) yl(2)],[0.5 0.5 0.5],'FaceAlpha',.25), hold on
text(-2.45,0.8*yl(2),'Baseline Period')
plot(time_fr(time_fr > -2.5),kldivs_inc)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')

subplot(2,3,5)
title("Dec: single-time FRs vs baseline")
plot(time_fr(time_fr > -2.5),kldivs_dec)
yl = get(gca,'YLim');
patch([-2.5 -1.5 -1.5 -2.5],[yl(1) yl(1) yl(2) yl(2)],[0.5 0.5 0.5],'FaceAlpha',.25), hold on
text(-2.45,0.8*yl(2),'Baseline Period')
plot(time_fr(time_fr > -2.5),kldivs_dec)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')


subplot(2,3,6)
title("Divergence b/t inc and dec neurons")
plot(time_fr(time_fr > -2.5),kldivs_cf)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')

saveas(gcf,[figures_path,'/KLdiv_avg'],'epsc')



%% KS Density Plots - full
figure('Position', [10 10 750 500])


% Force a minimum firing rate per-trial of 0.5
ix_inc = nanmean(incFRs_full(:,21:end),2) > 1;
ix_dec = nanmean(decFRs_full(:,21:end),2) > 1;


colors = [linspace(0.929,0.494,50);linspace(0.694,0.184,50);linspace(0.1250,0.556,50)]';

incFRs_full_base = incFRs_full(:,21:40);
decFRs_full_base = decFRs_full(:,21:40);
[f_base_inc,xi_base_inc] = ksdensity(incFRs_full_base(incFRs_full_base>0),0:0.5:50,'Support',[-0.001,250],'BoundaryCorrection','reflection');
[f_base_dec,xi_base_dec] = ksdensity(decFRs_full_base(decFRs_full_base>0),0:0.5:50,'Support',[-0.001,250],'BoundaryCorrection','reflection');

kldivs_inc = zeros(50,1);
kldivs_dec = zeros(50,1);
kldivs_cf = zeros(50,1);

colormap(colors)
for i=21:70
    subplot(2,3,1)
    [finc,xi] = ksdensity(incFRs_full(incFRs_full(:,i)>0,i),0:0.5:50,'Support',[-0.001,250],'BoundaryCorrection','reflection');
    kldivs_inc(i-20) = kldivergence(finc,f_base_inc);

    plot(xi,finc,'color',colors(i-20,:))
    hold on
    
    subplot(2,3,2)
    [fdec,xi]=ksdensity(decFRs_full(decFRs_full(:,i)>0,i),0:0.5:50,'Support',[-0.001,250],'BoundaryCorrection','reflection');
    kldivs_dec(i-20) = kldivergence(fdec,f_base_dec);
    plot(xi,fdec,'color',colors(i-20,:))
    hold on
    
    subplot(2,3,3)
    plot(xi,finc-fdec,'color',colors(i-20,:)), hold on
    kldivs_cf(i-20) = 0.5 * kldivergence(finc,fdec) + 0.5 * kldivergence(fdec,finc);
end


sgtitle('Density & Divergence of Individual Trials')
subplot(2,3,1)
title("Density of Inc Neuron FRs over time")
xlabel('Firing Rate'),xlim([0 40])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');
subplot(2,3,2)
title("Density of Dec Neuron FRs over time")
xlabel('Firing Rate'),xlim([0 20])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');
subplot(2,3,3)
title("Inc - Dec Neuron Densitities over time")
xlabel('Firing Rate'),xlim([0 20])
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'-2.5 s','-2 s','-1.5 s','-1 s','-0.5 s','0 s'},...
    'Direction','reverse');

subplot(2,3,4)
title("Inc: single-time Spike Counts vs baseline")
plot(time_fr(time_fr > -2.5),kldivs_inc)
yl = get(gca,'YLim');
patch([-2.5 -1.5 -1.5 -2.5],[yl(1) yl(1) yl(2) yl(2)],[0.5 0.5 0.5],'FaceAlpha',.25), hold on
text(-2.45,0.8*yl(2),'Baseline Period')
plot(time_fr(time_fr > -2.5),kldivs_inc)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')

subplot(2,3,5)
title("Dec: single-time Spike Counts vs baseline")
plot(time_fr(time_fr > -2.5),kldivs_dec)
yl = get(gca,'YLim');
patch([-2.5 -1.5 -1.5 -2.5],[yl(1) yl(1) yl(2) yl(2)],[0.5 0.5 0.5],'FaceAlpha',.25), hold on
text(-2.45,0.8*yl(2),'Baseline Period')
plot(time_fr(time_fr > -2.5),kldivs_dec)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')

subplot(2,3,6)
title("Divergence b/t inc and dec neurons")
plot(time_fr(time_fr > -2.5),kldivs_cf)
xlabel('Time'),xlim([-2.5,0])
ylabel('KL Divergence')


saveas(gcf,[figures_path,'/KLdiv_full'],'epsc')



%% Pie Charts for Frequency of Inc/Dec Neurons

%% Snap Main Plotting


%% Initial
clear all
close all

addpath(genpath('./utils'))

%% Set-up

nettype = 'flux'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['./data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);

% Loading for pre-allocation
load([fpath files(1).name],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks','time','time_fr')

decFRs_ga = nan(length(decFRs),length(files),3);
incFRs_ga = nan(length(incFRs),length(files),3);

incSpks_ga = nan(length(incSpks),length(files),3);
decSpks_ga = nan(length(decSpks),length(files),3);
nonSpks_ga = nan(length(nonSpks),length(files),3);

EEG_ga = nan(length(EEG),length(files),3);

[f,xi] = ksdensity(tx/1000,0:0.5:10);

tx_ga = nan(length(f),length(files),3);


fpath = ['./data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);


nInc = nan(length(files),3);
nDec = nan(length(files),3);
nNon = nan(length(files),3);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks','type')
    
    if sum(~isnan(tx)) < 15, continue, end
    
    incFRs_ga(:,ii,1) = incFRs - nanmean(incFRs(1:20));
    decFRs_ga(:,ii,1) = nanmean(decFRs - nanmean(decFRs(1:20,:),1),2);
    EEG_ga(:,ii,1) = EEG - nanmean(EEG(1:baseline));
    tx_ga(:,ii,1) = ksdensity(tx/1000,0:0.5:10);
    
    incSpks_ga(:,ii,1) = incSpks - nanmean(incSpks(1:baseline));
    decSpks_ga(:,ii,1) = decSpks - nanmean(decSpks(1:baseline));
    nonSpks_ga(:,ii,1) = nonSpks - nanmean(nonSpks(1:baseline));
    
    nInc(ii,1) = sum(strcmp(type,'inc'));
    nDec(ii,1) = sum(strcmp(type,'dec'));
    nNon(ii,1) = sum(strcmp(type,'non'));
    
    
end


% Add in the Ksyn 0.5
fpath = ['./data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks')
    
%     if all(isnan(tx)) || ,continue,end
    if sum(~isnan(tx)) < 10, continue, end
    
    incFRs_ga(:,ii,2) = incFRs - nanmean(incFRs(1:20));
    decFRs_ga(:,ii,2) = nanmean(decFRs - nanmean(decFRs(1:20,:),1),2);
    EEG_ga(:,ii,2) = EEG - nanmean(EEG(1:baseline));
    tx_ga(:,ii,2) = ksdensity(tx/1000,0:0.5:10);
    
    incSpks_ga(:,ii,2) = incSpks - nanmean(incSpks(1:baseline));
    decSpks_ga(:,ii,2) = decSpks - nanmean(decSpks(1:baseline));
    nonSpks_ga(:,ii,2) = nonSpks - nanmean(nonSpks(1:baseline));
    
    nInc(ii,2) = sum(strcmp(type,'inc'));
    nDec(ii,2) = sum(strcmp(type,'dec'));
    nNon(ii,2) = sum(strcmp(type,'non'));
    
end

% Add in the Ksyn 1
fpath = ['./data/' nettype '/Ksyn_' num2str(1) '/'];
files = dir([fpath,'*.mat']);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks')
    
    if sum(~isnan(tx)) < 15, continue, end
    
    incFRs_ga(:,ii,3) = incFRs - nanmean(incFRs(1:20));
    decFRs_ga(:,ii,3) = nanmean(decFRs - nanmean(decFRs(1:20,:),1),2);
    EEG_ga(:,ii,3) = EEG - nanmean(EEG(1:baseline));
    tx_ga(:,ii,3) = ksdensity(tx/1000,0:0.5:10);
    
    incSpks_ga(:,ii,3) = incSpks - nanmean(incSpks(1:baseline));
    decSpks_ga(:,ii,3) = decSpks - nanmean(decSpks(1:baseline));
    nonSpks_ga(:,ii,3) = nonSpks - nanmean(nonSpks(1:baseline));
    
    nInc(ii,3) = sum(strcmp(type,'inc'));
    nDec(ii,3) = sum(strcmp(type,'dec'));
    nNon(ii,3) = sum(strcmp(type,'non'));
    
end

nInc = nanmean(nInc,1);
nDec = nanmean(nDec,1);
nNon = nanmean(nNon,1);



%% Pie Charts for Frequency

figure

X = [nInc;nDec];



labels = {'Increasing','Decreasing'};

subplot(1,2,1)
p = pie([nInc,nDec,nNon]);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = 'red';
p(1).EdgeColor = 'white';

p(2).FaceColor = 'blue';
p(2).EdgeColor = 'white';

t(1).FontSize = 14;
t(2).FontSize = 14;


subplot(1,2,2)
p = pie([nInc,nDec]);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = 'red';
p(1).EdgeColor = 'white';

p(2).FaceColor = 'blue';
p(2).EdgeColor = 'white';

t(1).FontSize = 14;
t(2).FontSize = 14;

%% Pie Chart taking into account ramping behavior


figure

X = [nInc;nDec];


% nInc_step = sum(incGain > nanmean(incGain));
% nDec_step = sum(decGain > nanmean(decGain));
% 
% nInc_ramp = sum(incGain < nanmean(incGain));
% nDec_ramp = sum(decGain < nanmean(decGain));

nInc_step = sum(incGain > 0.75);
nDec_step = sum(decGain > 0.75);

nInc_ramp = sum(incGain < 0.75);
nDec_ramp = sum(decGain < 0.75);



labels = {'Increasing','Decreasing'};

subplot(1,2,1)
p = pie([nInc_step,nDec_step]);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = 'red';
p(1).EdgeColor = 'white';

p(2).FaceColor = 'blue';
p(2).EdgeColor = 'white';

t(1).FontSize = 14;
t(2).FontSize = 14;


subplot(1,2,2)
p = pie([nInc_ramp,nDec_ramp]);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = 'red';
p(1).EdgeColor = 'white';

p(2).FaceColor = 'blue';
p(2).EdgeColor = 'white';

t(1).FontSize = 14;
t(2).FontSize = 14;

%% Violin plots for taus

tau_cutoff = 100;
nInc2 = length(incTaus(incTaus < tau_cutoff));
nDec2 = length(decTaus(decTaus < tau_cutoff));
taumat = nan(max([nInc2,nDec2]),2);



taumat(1:nInc2,1) = incTaus(incTaus < tau_cutoff);
taumat(1:nDec2,2) = decTaus(decTaus < tau_cutoff);

violin(taumat,'xlabel',{'inc','dec','nochange'},...
    'facecolor',[1 0 0;0 0 1;0 0 0]);



figure('position',[0 0 200 200])

histogram(incTaus(incTaus < tau_cutoff),...
    'FaceColor','red','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:5:tau_cutoff,'Normalization','probability')
hold on
histogram(decTaus(decTaus < tau_cutoff),...
    'FaceColor','blue','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:5:tau_cutoff,'Normalization','probability')
xticks([0 50 100])
xlabel('Autocorrelation \tau','FontSize',16)
ylabel('Proportion','FontSize',16)

saveas(gcf,['figures/real_histTaus.svg'])

%% Separate Y-Axes for increasing & decreasing neurons

figure('Position',[0 0 400 320])
hold on

ymax = max(incSpks,[],'all');

gain_cutoff = 0.4;

% Increasing neurons
% p1K0 = plot_signal_ci(time,incSpks(incGain <= nanmean(incGain),:),'Color',[1 0.3 0.2]); % Ksyn = 0
% p1K5 = plot_signal_ci(time,incSpks(incGain >= nanmean(incGain),:),'Color',[1 0.3 0.2]); % Ksyn = 0.5

p1K0 = plot_signal_ci(time,incSpks(incGain <= gain_cutoff & ismember(incRegion,mfcRegions),:),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(time,incSpks(incGain >= gain_cutoff & ismember(incRegion,mfcRegions),:),'Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

yyaxis left

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis(1).FontSize = 14;
ax.YAxis(2).FontSize = 14;

ylabel('\DeltaFiring Rate (Hz)','FontSize',20)


% Setting ylimit so they work together
yl = ylim;
ylim([yl(1) - 2, yl(2) + 1])
% % if strcmp(nettype,'flux')
%     ylim([yl(1) - 15, yl(2)])
% else
%     ylim([yl(1) - 50, yl(2)])
% end
yl = ylim;

set(gca,'ycolor',[1 0.3 0.2])

% decreasing neurons

yyaxis right
% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% p2K0 = plot_signal_ci(time,decSpks(decGain <= nanmean(decGain),:),'Color',[0 0.5 1]); % Ksyn = 0
% p2K5 = plot_signal_ci(time,decSpks(decGain >= nanmean(decGain),:),'Color',[0 0.5 1]); % Ksyn = 0.5

p2K0 = plot_signal_ci(time,decSpks(decGain <= gain_cutoff & ismember(decRegion,mfcRegions),:),'Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(time,decSpks(decGain >= gain_cutoff & ismember(decRegion,mfcRegions),:),'Color',[0 0.5 1]); % Ksyn = 0.5
set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h.Position(1) = 0.5;
% label_h.Position(2) = -0.;

ylim(yl/2)
% if strcmp(nettype,'flux')
%     
% else
%     ylim(yl/4)
% end

% ax.YColor = [0 0.5 1];
set(gca,'ycolor',[0 0.5 1])


yline(0,'--k','LineWidth',2)


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Inc, Steplike','Inc, Slow Ramping',...
    'Dec, Steplike','Dec, Slow Ramping')

set(legend,'Location','NorthWest','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time rel. movement (s)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)
xlim([-3 0])





% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/real_FRs_both_',num2str(gain_cutoff),'_mfc.svg'])
saveas(gcf,['figures/real_FRs_both_',num2str(gain_cutoff),'_mfc.png'])


%% Autocorrelation

figure('Position',[0 0 350 320])
hold on


lags = 0:1000;

ACFsI = incACF;
ACFsD = decACF;

ACFsI(:,lags <= 10,:) = nan;
ACFsD(:,lags <= 10,:) = nan;

lags(lags <=10) = nan;

ACFsI = ACFsI(:,~isnan(lags),:);
ACFsD = ACFsD(:,~isnan(lags),:);

lags = lags(~isnan(lags));

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
yline(0,'--k','LineWidth',2)

% Increasing neurons
p1K0 = plot_signal_ci(lags,ACFsI(incGain >= gain_cutoff & ismember(incRegion,mfcRegions),:),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(lags,ACFsD(decGain >= gain_cutoff & ismember(decRegion,mfcRegions),:),'Color',[0 0.5 1]); % Ksyn = 0.5
% set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
xlim([0 1000])
ylim([-0.01 0.055])
% xlim([0 500])
% set(gca,'ycolor','r')


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5],...
    'Increasing','Decreasing')

set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/real_ACFs_mfc.svg'])
saveas(gcf,['figures/real_ACFs_mfc.png'])
%%
figure('Position',[0 0 270 250])
hold on


lags = 0:1000;

ACFsI = incACF;
ACFsD = decACF;

ACFsI(:,lags <= 10,:) = nan;
ACFsD(:,lags <= 10,:) = nan;

lags(lags <=10) = nan;

ACFsI = ACFsI(:,~isnan(lags),:);
ACFsD = ACFsD(:,~isnan(lags),:);

lags = lags(~isnan(lags));

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
yline(0,'--k','LineWidth',2)

% Increasing neurons
p2K0 = plot_signal_ci(lags,ACFsI(incGain > nanmedian(incGain),:),'Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(lags,ACFsI(incGain <= nanmedian(incGain),:),'Color',[0 0.5 1]); % Ksyn = 0.5
set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')
% Increasing neurons
p1K0 = plot_signal_ci(lags,ACFsD(decGain > nanmedian(decGain),:),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(lags,ACFsD(decGain <= nanmedian(decGain),:),'Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K5,'LineStyle','--')

ylabel('Autocorrelation','FontSize',20)
xlim([0 1000])
ylim([-0.005 0.035])
% xlim([0 500])
% set(gca,'ycolor','r')


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5],...
    'Inc, High Gain','Inc, Low Gain', 'Dec, High Gain','Dec, Low Gain')

set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/real_ACFs_gainsplit.svg'])
saveas(gcf,['figures/real_ACFs_gainsplit.png'])


%% Calculating onset of ramping

% Idea: find first time at which rank-sum test between data & baseline is
% significant.
gain_cutoff=0.4;
slowrampix = find(incGain >= gain_cutoff & ismember(incRegion,mfcRegions));
time2use = (time > -3.5 & time < -1);
incOnsets = zeros(length(slowrampix),1);
for ii = 1:length(slowrampix)
%     sigs = find(incSpks_sig(slowrampix(ii),:,1) < 0.01 & incSpks_sig(slowrampix(ii),:,2) > 0);
%     incOnsets(ii) = time(sigs(1));
    incOnsets(ii) = time(findchangepts(incSpks(slowrampix(ii),time2use),'Statistic','linear'));
end

slowrampix = find(incGain >= gain_cutoff & ismember(incRegion,mfcRegions));

decOnsets = zeros(length(slowrampix),1);
for ii = 1:length(slowrampix)
%     sigs = find(decSpks_sig(slowrampix(ii),:,1) < 0.01 & decSpks_sig(slowrampix(ii),:,2) < 0);
%     decOnsets(ii) = time(sigs(1));
    decOnsets(ii) = time(findchangepts(decSpks(slowrampix(ii),time2use),'Statistic','linear'));
end

figure('Position',[0 0 400 300])
histogram(incOnsets,...
    'FaceColor',[1 0.3 0.2],'LineWidth',1.5,'EdgeAlpha',1,'FaceAlpha',0.5,...
    'BinEdges',-3.5:0.25:0)
hold on
histogram(decOnsets,...
    'FaceColor',[0 0.5 1],'LineWidth',1.5,'EdgeAlpha',1,'FaceAlpha',0.5,...
    'BinEdges',-3.5:0.25:0)

ax = gca;
axis off
ax.YAxis.Visible = 'on';
ax.XAxis.Visible = 'on';
ax.FontSize = 16; 
ylabel('Frequency','FontSize',20) 
xlabel('Ramping Onset (s)','FontSize',20) 


saveas(gcf,['figures/ramponset_mfc.svg'])
saveas(gcf,['figures/ramponset_mfc.png'])



% incOnsets = zeros(size(incSpks,1),1);
% for ii = 1:nInc
% %     sigs = find(incSpks_sig(slowrampix(ii),:,1) < 0.01 & incSpks_sig(slowrampix(ii),:,2) > 0);
% %     incOnsets(ii) = time(sigs(1));
%     incOnsets(ii) = time(findchangepts(incSpks(ii,time2use),'Statistic','linear'));
% end
% 
% decOnsets = zeros(size(decSpks,1),1);
% for ii = 1:nDec
% %     sigs = find(decSpks_sig(slowrampix(ii),:,1) < 0.01 & decSpks_sig(slowrampix(ii),:,2) < 0);
% %     decOnsets(ii) = time(sigs(1));
%     decOnsets(ii) = time(findchangepts(decSpks(ii,time2use),'Statistic','linear'));
% end


% histogram(incOnsets(ismember(incRegion,mfcRegions)),...
%     'FaceColor','red','EdgeAlpha',1,'FaceAlpha',0.5,...
%     'BinEdges',-3.5:0.25:0)
% hold on
% histogram(decOnsets(ismember(decRegion,mfcRegions)),...
%     'FaceColor','blue','EdgeAlpha',1,'FaceAlpha',0.5,...
%     'BinEdges',-3.5:0.25:0)


