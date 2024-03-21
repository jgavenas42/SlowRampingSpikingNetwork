%% Snap Main Plotting


%% Initial
clear all
close all

addpath(genpath('./utils'))

%% Set-up

nettype = 'flux'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['../data/Ksyn_cf/' nettype '/']; % 4 clusters
% fpath = ['../data/2clusters/' nettype '/']; % 2, 3, or 5 clusters

files = dir([fpath,'*.mat']);

% Loading for pre-allocation
load([fpath files(1).name],'incSpks','decSpks','time','Ksyn',...
    'EEG','acfs_all','acfs_clust')

decSpks_ga = nan(length(time),length(files));
incSpks_ga = nan(length(time),length(files));

eeg_ga = nan(length(time),length(files));

acfsN_ga = nan(1001,length(files));
acfsC_ga = nan(1001,length(files));

Ksyns_ga = nan(length(files),1);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'incSpks','decSpks','Ksyn','EEG','acfs_all','acfs_clust')
    
    
    incSpks_ga(:,ii) = incSpks - nanmean(incSpks(1:baseline));
    decSpks_ga(:,ii) = decSpks - nanmean(decSpks(1:baseline));
    
    EEG = smoothdata(EEG,2,'gaussian',200);
    
    eeg_ga(:,ii) = nanmean(EEG - nanmean(EEG(:,1:baseline),2),1);
    
    acfsN_ga(:,ii) = nanmean(acfs_all,1);
    acfsC_ga(:,ii) = nanmean(acfs_clust,1);
    
    Ksyns_ga(ii) = Ksyn;
end
    

%% Average across networks with the same Ksyn

Ksyns = 0:0.05:0.5;

incSpks_all = nan(length(time),length(Ksyns));
decSpks_all = nan(length(time),length(Ksyns));

acfsN = nan(1001,length(Ksyns));
acfsC = nan(1001,length(Ksyns));

eeg = nan(length(time),length(Ksyns));

for ii = 1:length(Ksyns)
    
    incSpks_all(:,ii) = nanmean(incSpks_ga(:,Ksyns_ga == Ksyns(ii)),2); 
    decSpks_all(:,ii) = nanmean(decSpks_ga(:,Ksyns_ga == Ksyns(ii)),2); 
    
    eeg(:,ii) = nanmean(eeg_ga(:,Ksyns_ga == Ksyns(ii)),2);
    
    acfsN(:,ii) = nanmean(acfsN_ga(:,Ksyns_ga == Ksyns(ii)),2);
    acfsC(:,ii) = nanmean(acfsC_ga(:,Ksyns_ga == Ksyns(ii)),2);
    
end

%% Plot with varying Ksyn




figure('Position',[0 0 400 320])
hold on
% colormap(spring(length(Ksyns)));
% cmap = spring(length(Ksyns));
cmap = brewermap(ceil(length(Ksyns)*1.2),'-YlOrRd');

for ii = 1:length(Ksyns)
    plot(time/1000,incSpks_all(:,ii),'Color',cmap(ii,:),'LineWidth',2)
end

xlim([-3 0])
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('\Delta Firing Rates (Hz)','FontSize',20)
% title('Increasing Neurons','FontSize',24)

colormap(cmap(1:11,:));
a = colorbar;
% caxis([0 0.])

a.Ticks = [0 0.25 0.5];
a.TickLabels = {'0%','25%','50%'};
a.FontSize = 14;

ylabel(a,'% Slow Synapses','FontSize',20,'Rotation',270);
a.Label.Position(1) = 5;
a.Label.Position(2) = 0.25;



caxis([0 0.5])


% saveas(gcf,['figures/' nettype '_KSynCF_inc.png'])
% saveas(gcf,['figures/' nettype '_KSynCF_inc.svg'])


figure('Position',[0 0 400 320])

hold on
% colormap(flip(cool(length(Ksyns))));
% cmap = flip(cool(length(Ksyns)));
% cmap = winter(length(Ksyns));
cmap = brewermap(ceil(length(Ksyns)*1.2),'-PuBuGn');
for ii = 1:length(Ksyns)
    plot(time/1000,decSpks_all(:,ii),'Color',cmap(ii,:),'LineWidth',2)
end
xlim([-3 0])

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('\Delta Firing Rates (Hz)','FontSize',20)
% title('Decreasing Neurons','FontSize',24)
colormap(cmap(1:11,:));


b = colorbar;


b.Ticks = [0 0.25 0.5];
b.TickLabels = {'0%','25%','50%'};
b.FontSize = 14;

ylabel(b,'% Slow Synapses','FontSize',20,'Rotation',270);
b.Label.Position(1) = 5;
b.Label.Position(2) = 0.25;

caxis([0 0.5])
% saveas(gcf,['figures/' nettype '_KSynCF_dec.svg'])
% saveas(gcf,['figures/' nettype '_KSynCF_dec.png'])


%% EEG



figure('Position',[0 0 400 320])
hold on

% cmap = flip(brewermap(length(Ksyns)*2,'Greens'));
% cmap = cmap(1:11,:);

% cmap = parula(length(Ksyns)*1.5);
cmap = brewermap(ceil(length(Ksyns)*1.5),'-YlGn');

for ii = 1:length(Ksyns)
    plot(time/1000,eeg(:,ii),'Color',cmap(ii,:),'LineWidth',2)
end

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('Voltage (mV)','FontSize',20)
% xlim([-2 0.5])
% title('Increasing Neurons','FontSize',24)
colormap(cmap(1:11,:));
a = colorbar;
a.Ticks = [0 0.25 0.5];
a.TickLabels = {'0%','25%','50%'};
a.FontSize = 14;

ylabel(a,'% Slow Synapses','FontSize',20,'Rotation',270);
a.Label.Position(1) = 5;
a.Label.Position(2) = 0.25;

caxis([0 0.5])
xlim([-3 0])
% ylim([-0.03,0.005])

% saveas(gcf,['figures/' nettype '_KSynCF_eeg.svg'])
% saveas(gcf,['figures/' nettype '_KSynCF_eeg.png'])


%% ACF Neurons

acfsT = flip(acfsN);

lags = -1000:1000;

acfsT = cat(1,acfsT(1:1000,:),acfsN);

acfsT(lags >= -10 & lags <= 10,:) = nan;
lags(lags >= -10 & lags <= 10) = nan;


figure('Position',[0 0 500 400])
hold on
% cmap = parula(length(Ksyns));

% cmap = [linspace(0,1,length(Ksyns));zeros(length(Ksyns),1)';linspace(0,1,length(Ksyns))]';

cmap = flip(brewermap(length(Ksyns)*2,'PuBu'));

% cmap = autumn(length(Ksyns)*2);

cmap = cmap(1:11,:);

for ii = 1:length(Ksyns)
    plot(lags,acfsT(:,ii),'Color',cmap(ii,:),'LineWidth',2)
end
xlabel('Time Lag (ms)','FontSize',20)
ylabel('Autocorrelation','FontSize',20)
xlim([-1000 1000])
ylim([-0.01 0.04])
% title('Increasing Neurons','FontSize',24)
colormap(cmap);
a = colorbar;
ylabel(a,'Ksyn','FontSize',20,'Rotation',270);
a.Label.Position(1) = 4;
a.Label.Position(2) = 0.25;

caxis([0 0.5])

% saveas(gcf,['figures/' nettype '_KSynCF_acfsN.svg'])
% saveas(gcf,['figures/' nettype '_KSynCF_acfsN.png'])

%% ACF Clusters

acfsT = flip(acfsC);

lags = -1000:1000;

acfsT = cat(1,acfsT(1:1000,:),acfsC);

acfsT(lags >= -10 & lags <= 10,:) = nan;
lags(lags >= -10 & lags <= 10) = nan;


figure('Position',[0 0 500 400])
hold on
% cmap = parula(length(Ksyns));

% cmap = [linspace(0,1,length(Ksyns));zeros(length(Ksyns),1)';linspace(0,1,length(Ksyns))]';

cmap = flip(brewermap(length(Ksyns)*2,'RdPu'));

cmap = cmap(1:11,:);

for ii = 1:length(Ksyns)
    plot(lags,acfsT(:,ii),'Color',cmap(ii,:),'LineWidth',2)
end

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

xlabel('Time Lag (ms)','FontSize',20)
ylabel('Autocorrelation','FontSize',20)
% xlim([0 300])
% ylim([-0.06 0.04])
% title('Increasing Neurons','FontSize',24)
colormap(cmap);
a = colorbar;

a.Ticks = [0 0.25 0.5];
a.TickLabels = {'0%','25%','50%'};
a.FontSize = 14;

ylabel(a,'% Slow Synapses','FontSize',20,'Rotation',270);
a.Label.Position(1) = 5;
a.Label.Position(2) = 0.25;



caxis([0 0.5])

% saveas(gcf,['figures/' nettype '_KSynCF_acfsC.svg'])
% saveas(gcf,['figures/' nettype '_KSynCF_acfsC.png'])


    
