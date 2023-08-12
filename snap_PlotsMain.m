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
    'incSpks','decSpks','nonSpks','time','time_fr','taus')

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

incTaus = nan(320,length(files),3);
decTaus = nan(320,length(files),3);

incACFs = nan(length(files),1001,3);
decACFs = nan(length(files),1001,3);

incACFs_clust = nan(length(files),1001,3);
decACFs_clust = nan(length(files),1001,3);

ACFs_neuron = nan(length(files),1001,3);
ACFs_clust = nan(length(files),1001,3);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks','type',...
    'acfs_all','acfs_clust','max_clust')
    
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
    
    incTaus(1:nInc(ii,1),ii,1) = taus(strcmp(type,'inc'));
    decTaus(1:nDec(ii,1),ii,1) = taus(strcmp(type,'dec'));
    
    incACFs(ii,:,1) = nanmean(acfs_all(strcmp(type,'inc'),:),1);
    decACFs(ii,:,1) = nanmean(acfs_all(strcmp(type,'dec'),:),1);
    
    incACFs_clust(ii,:,1) = acfs_clust(max_clust,:);
    decACFs_clust(ii,:,1) = nanmean(acfs_clust(setdiff(1:4,max_clust),:),1);
    
    ACFs_neuron(ii,:,1) = nanmean(acfs_all,1);
    ACFs_clust(ii,:,1) = nanmean(acfs_clust,1);

    

    
    
    
    
end


% Add in the Ksyn 0.5
fpath = ['./data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks','taus',...
    'acfs_all','acfs_clust','max_clust')
    
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
    
    incTaus(1:nInc(ii,2),ii,2) = taus(strcmp(type,'inc'));
    decTaus(1:nDec(ii,2),ii,2) = taus(strcmp(type,'dec'));
    
    incACFs(ii,:,2) = nanmean(acfs_all(strcmp(type,'inc'),:),1);
    decACFs(ii,:,2) = nanmean(acfs_all(strcmp(type,'dec'),:),1);
    
    incACFs_clust(ii,:,2) = acfs_clust(max_clust,:);
    decACFs_clust(ii,:,2) = nanmean(acfs_clust(setdiff(1:4,max_clust),:),1);
    
    ACFs_neuron(ii,:,2) = nanmean(acfs_all,1);
    ACFs_clust(ii,:,2) = nanmean(acfs_clust,1);
    
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
    
    incTaus(1:nInc(ii,3),ii,3) = taus(strcmp(type,'inc'));
    decTaus(1:nDec(ii,3),ii,3) = taus(strcmp(type,'dec'));
    
end

nInc = nanmean(nInc,1);
nDec = nanmean(nDec,1);
nNon = nanmean(nNon,1);

%% Increasing & Decreasing Neuron FRs
figure
hold on

% Plot the Ksyn = 0
p1K0 = plot_signal_ci(time/1000,incSpks_ga(:,:,1)','r-');
p2K0 = plot_signal_ci(time/1000,decSpks_ga(:,:,1)','b-');

% Plot Ksyn = 0.5
p1K5 = plot_signal_ci(time/1000,incSpks_ga(:,:,2)','r--');
p2K5 = plot_signal_ci(time/1000,decSpks_ga(:,:,2)','b--');
set(p1K5,'LineStyle','--')
set(p2K5,'LineStyle','--')


% Plot Ksyn = 1
p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:');
p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:');
set(p1K1,'LineStyle',':')
set(p2K1,'LineStyle',':')

% plot_signal_ci(time/1000,nonSpks_ga','k')
legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
    'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
    'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')

set(legend,'Location','NorthWest','FontSize',16)
legend('boxoff')

grid
xlabel('Time (s); 0 = TX','FontSize',20)
ylabel('Baselined Firing Rate (Hz)','FontSize',20)
title('Early Ramping in Slow Synapse Nets','FontSize',24)

if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
% saveas(gcf,['figures/' nettype '_FRs.svg'])

%% Separate Y-Axes for increasing & decreasing neurons

figure('Position',[0 0 400 320])
hold on

ymax = max(incSpks_ga,[],'all');

% Increasing neurons
p1K0 = plot_signal_ci(time/1000,incSpks_ga(:,:,1)','Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(time/1000,incSpks_ga(:,:,2)','Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

yyaxis left

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis(1).FontSize = 14;
ax.YAxis(2).FontSize = 14;

ylabel('\DeltaFiring Rate (Hz)','FontSize',20)
set(gca,'ycolor',[1 0.3 0.2])

% Setting ylimit so they work together
yl = ylim;

if strcmp(nettype,'flux')
    ylim([yl(1) - 15, yl(2)])
else
    ylim([yl(1) - 50, yl(2)])
end
yl = ylim;


% decreasing neurons

yyaxis right
p2K0 = plot_signal_ci(time/1000,decSpks_ga(:,:,1)','Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(time/1000,decSpks_ga(:,:,2)','Color',[0 0.5 1]); % Ksyn = 0.5
set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
set(gca,'ycolor',[0 0.5 1])

if strcmp(nettype,'flux')
    ylim(yl/4)
else
    ylim(yl/4)
end
yline(0,'--k','LineWidth',2)


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Inc, 100% Fast','Inc, 50% Slow',...
    'Dec, 100% Fast','Dec, 50% Slow')

set(legend,'Location','NorthWest','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time rel. threshold-crossing (s)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)
xlim([-3 0])


% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_FRs_both.svg'])
saveas(gcf,['figures/' nettype '_FRs_both.png'])


%% Increasing & Decreasing FRs Separately
figure
hold on

% Plot the Ksyn = 0
p1K0 = plot_signal_ci(time/1000,incSpks_ga(:,:,1)','r');
set(p1K0,'Color',[1 0.3 0.2])

% Plot Ksyn = 0.5
p1K5 = plot_signal_ci(time/1000,incSpks_ga(:,:,2)','r');
set(p1K5,'LineStyle','--','Color',[1 0.3 0.2])


% Plot Ksyn = 1
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r');
% set(p1K1,'LineStyle',':')

% plot_signal_ci(time/1000,nonSpks_ga','k')
legend([p1K0 p1K5 p1K1 ],...
    'Ksyn=0.0','Ksyn=0.5','Ksyn=1')

set(legend,'Location','NorthWest','FontSize',16)
legend('boxoff')

grid
xlabel('Time (s); 0 = TX','FontSize',20)
ylabel('Baselined Firing Rate (Hz)','FontSize',20)
title('Early Ramping in Slow Synapse Nets','FontSize',24)


% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_FRs_inc.svg'])


figure
hold on

% Plot the Ksyn = 0
p2K0 = plot_signal_ci(time/1000,decSpks_ga(:,:,1)','b');

% Plot Ksyn = 0.5
p2K5 = plot_signal_ci(time/1000,decSpks_ga(:,:,2)','b');
set(p2K5,'LineStyle','--')


% Plot Ksyn = 1
p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b');
set(p2K1,'LineStyle',':')

% plot_signal_ci(time/1000,nonSpks_ga','k')
legend([p2K0 p2K5  p2K1],...
    'Ksyn=0.0','Ksyn=0.5','Ksyn=1')

set(legend,'Location','SouthWest','FontSize',16)
legend('boxoff')

grid
xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('Baselined Firing Rate (Hz)','FontSize',20)
title('Inhibitory Neurons','FontSize',24)


% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
% saveas(gcf,['figures/' nettype '_FRs_dec.svg'])
%% EEG

figure('Position',[0 0 400 320])

hold on
p1 = plot_signal_ci(time/1000,EEG_ga(:,:,1)','Color',[0 0.5 0]);
p2 = plot_signal_ci(time/1000,EEG_ga(:,:,2)','Color',[0 0.5 0]);
% p3 = plot_signal_ci(time/1000,EEG_ga(:,:,3)','g');

set(p1,'LineStyle','-')
set(p2,'LineStyle','--')
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlim([-3 0.])
% set(p3,'LineStyle',':','Color',[0 0.5 0])

% legend([p1 p2 p3],...
%     'Ksyn=0.0','Ksyn=0.5','Ksyn=1')
legend([p1 p2],...
    '100% Fast','50% Slow')

set(legend,'Location','SouthWest','FontSize',16)
legend('boxoff')


grid off
xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('Voltage (mV)','FontSize',20)
ylim([-0.03,0.005])
% title('RP in Slow Synapse Nets','FontSize',24)
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_EEG.svg'])
saveas(gcf,['figures/' nettype '_EEG.png'])


%% TX Distribution

figure
hold on
p1 = plot_signal_ci(xi,tx_ga(:,:,1)','Color','black');
p2 = plot_signal_ci(xi,tx_ga(:,:,2)','Color','black');
set(p2,'LineStyle','--')
% p3 = plot_signal_ci(xi,tx_ga(:,:,3)','k:');

legend([p1 p2],...
    'Ksyn=0.0','Ksyn=0.5')
set(legend,'Location','NorthEast','FontSize',16)
legend('boxoff')


grid off
xlabel('Threshold-Crossing Times','FontSize',20)
ylabel('Density','FontSize',20)
% title('Threshold-Crossing Distribution','FontSize',24)
saveas(gcf,['figures/' nettype '_TX.svg'])
saveas(gcf,['figures/' nettype '_TX.png'])



%% Pie Charts for Frequency

figure

X = [nInc;nDec];
explode = [2 0];

edgecolor = 'black';
edgewidth = 2;



labels = {'Increasing','Decreasing'};

subplot(1,3,1)
p = pie(X(:,1),explode);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = [1 0.3 0.2];
p(1).EdgeColor = edgecolor;
p(1).LineWidth = edgewidth;


p(2).FaceColor = [0 0.5 1];
p(2).EdgeColor = edgecolor;
p(2).LineWidth = edgewidth;

t(1).FontSize = 14;
t(2).FontSize = 14;
% title('Ksyn = 0.0','FontSize',20)


subplot(1,3,2)
p = pie(X(:,2),explode);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = [1 0.3 0.2];
p(1).EdgeColor = edgecolor;
p(1).LineWidth = edgewidth;

p(2).FaceColor = [0 0.5 1];
p(2).EdgeColor = edgecolor;
p(2).LineWidth = edgewidth;


t(1).FontSize = 14;
t(2).FontSize = 14;

% title('Ksyn = 0.5','FontSize',20)


subplot(1,3,3)
% Take the numbers from the real data
% nInc = 67, nDec = 181
% restricting to gain < 0.6, 54 and 130
% p = pie([67,181],explode);
p = pie([54, 130],explode);
t = p(2:2:end);
p = p(1:2:end);

p(1).FaceColor = [1 0.3 0.2];
p(1).EdgeColor = edgecolor;
p(1).LineWidth = edgewidth;

p(2).FaceColor = [0 0.5 1];
p(2).EdgeColor = edgecolor;
p(2).LineWidth = edgewidth;


t(1).FontSize = 14;
t(2).FontSize = 14;

% title('Real Data','FontSize',20)

% legend(labels)

saveas(gcf,['figures/' nettype '_typeFreqs.svg'])


%% Taus Comparison

ksyn2use = 0.5;

if ksyn2use == 0
    ix2use = 1;
elseif ksyn2use == 0.5
    ix2use = 2;
elseif ksyn2use == 1
    ix2use = 3;
else
    error('Choose a better ksyn')
end

incTaus_this = incTaus(:,:,ix2use);
decTaus_this = decTaus(:,:,ix2use);

tau_cutoff = 100;
% nInc2 = length(incTaus(incTaus_this < tau_cutoff));
% nDec2 = length(decTaus(decTaus_this < tau_cutoff));
% taumat = nan(max([nInc2,nDec2]),2);
% 
% 
% 
% taumat(1:nInc2,1) = incTaus_this(incTaus_this < tau_cutoff);
% taumat(1:nDec2,2) = decTaus_this(decTaus_this < tau_cutoff);
% 
% violin(taumat,'xlabel',{'inc','dec','nochange'},...
%     'facecolor',[1 0 0;0 0 1;0 0 0]);

figure('position',[0 0 200 200])

histogram(incTaus_this(incTaus_this < tau_cutoff),...
    'FaceColor','red','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:5:tau_cutoff,'Normalization','probability')
hold on
histogram(decTaus_this(decTaus_this < tau_cutoff),...
    'FaceColor','blue','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:5:tau_cutoff,'Normalization','probability')
xticks([0 50 100])
xlabel('Autocorrelation \tau','FontSize',16)
ylabel('Proportion','FontSize',16)


saveas(gcf,['figures/' nettype '_' num2str(ksyn2use) '_histTaus.svg'])


%% Autocorrelation Neurons

figure('Position',[0 0 270 250])
hold on

incACFs_flip = flip(incACFs,2);
decACFs_flip = flip(decACFs,2);

incACFs_t = cat(2,incACFs_flip(:,1:1000,:),incACFs);
decACFs_t = cat(2,decACFs_flip(:,1:1000,:),decACFs);

lags = -1000:1000;

incACFs_t(:,lags >= -10 & lags <= 10,:) = nan;
decACFs_t(:,lags >= -10 & lags <= 10,:) = nan;

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
yline(0,'--k','LineWidth',2)
% decreasing neurons

% yyaxis right
p2K0 = plot_signal_ci(lags,decACFs_t(:,:,1),'Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(lags,decACFs_t(:,:,2),'Color',[0 0.5 1]); % Ksyn = 0.5
set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')
% Increasing neurons
p1K0 = plot_signal_ci(lags,incACFs_t(:,:,1),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(lags,incACFs_t(:,:,2),'Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
ylim([-0.02 0.08])
xlim([0 1000])


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Inc, 100% Fast','Inc, 50% Slow',...
    'Dec, 100% Fast','Dec, 50% Slow')

set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_ACFs_neurons.svg'])
saveas(gcf,['figures/' nettype '_ACFs_neurons.png'])

%% Autocorrelation Neurons slow ramping only
figure('Position',[0 0 350 320])
hold on

incACFs_flip = flip(incACFs,2);
decACFs_flip = flip(decACFs,2);

incACFs_t = cat(2,incACFs_flip(:,1:1000,:),incACFs);
decACFs_t = cat(2,decACFs_flip(:,1:1000,:),decACFs);

lags = -1000:1000;

incACFs_t(:,lags >= -10 & lags <= 10,:) = nan;
decACFs_t(:,lags >= -10 & lags <= 10,:) = nan;

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
yline(0,'--k','LineWidth',2)
% decreasing neurons

% yyaxis right
% p2K0 = plot_signal_ci(lags,decACFs_t(:,:,1),'Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(lags,decACFs_t(:,:,2),'Color',[0 0.5 1]); % Ksyn = 0.5
% set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')
% Increasing neurons
% p1K0 = plot_signal_ci(lags,incACFs_t(:,:,1),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(lags,incACFs_t(:,:,2),'Color',[1 0.3 0.2]); % Ksyn = 0.5
% set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
ylim([-0.005 0.025])
xlim([0 1000])


% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([ p1K5  p2K5 ],...
    'Increasing',...
    'Deccreasing')

set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_ACFs_neurons_slowramponly.svg'])
saveas(gcf,['figures/' nettype '_ACFs_neurons_slowramponly.png'])


%% Autocorrelation Clusters
% 
% figure('Position',[0 0 450 400])
% hold on
% 
% incACFs_flip = flip(incACFs_clust,2);
% decACFs_flip = flip(decACFs_clust,2);
% 
% incACFs_t = cat(2,incACFs_flip(:,1:1000,:),incACFs_clust);
% decACFs_t = cat(2,decACFs_flip(:,1:1000,:),decACFs_clust);
% 
% lags = -1000:1000;
% 
% incACFs_t(:,lags >= -10 & lags <= 10,:) = nan;
% decACFs_t(:,lags >= -10 & lags <= 10,:) = nan;
% 
% 
% % Increasing neurons
% p1K0 = plot_signal_ci(lags,incACFs_t(:,:,1),'r-'); % Ksyn = 0
% p1K5 = plot_signal_ci(lags,incACFs_t(:,:,2),'r--'); % Ksyn = 0.5
% set(p1K5,'LineStyle','--')
% % p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% % set(p1K1,'LineStyle',':')
% 
% ylabel('Autocorrelation','FontSize',20)
% % ylim([-0.1 0.2])
% % xlim([0 500])
% % set(gca,'ycolor','r')
% 
% 
% 
% 
% % decreasing neurons
% 
% % yyaxis right
% p2K0 = plot_signal_ci(lags,decACFs_t(:,:,1),'b-'); % Ksyn = 0
% p2K5 = plot_signal_ci(lags,decACFs_t(:,:,2),'b--'); % Ksyn = 0.5
% set(p2K5,'LineStyle','--')
% % p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% % set(p2K1,'LineStyle',':')
% 
% % label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% % label_h.Position(1) = 0.9;
% % label_h.Position(2) = 1;
% % set(gca,'ycolor','b')
% 
% 
% 
% % legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
% %     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
% %     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
% legend([p1K0 p1K5 p2K0 p2K5 ],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5')
% 
% set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
% legend('boxoff')
% 
% grid off
% xlabel('Time Lag (ms)','FontSize',20)
% % title('Early Ramping in Slow Synapse Nets','FontSize',24)
% 
% 
% 
% % if strcmp(nettype,'flux'),ylim([-6 26]),end
% % plot_signal_ci(time/1000,EEG_ga' * 500,'g')
% saveas(gcf,['figures/' nettype '_ACFs_clusters.svg'])
% saveas(gcf,['figures/' nettype '_ACFs_clusters.png'])
% 
% 
% %% Autocorrelation Neurons
% 
% figure('Position',[0 0 450 400])
% hold on
% 
% 
% lags = 0:1000;
% 
% 
% % Increasing neurons
% p1K0 = plot_signal_ci(lags,incACFs(:,:,1),'r-'); % Ksyn = 0
% p1K5 = plot_signal_ci(lags,incACFs(:,:,2),'r--'); % Ksyn = 0.5
% set(p1K5,'LineStyle','--')
% % p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% % set(p1K1,'LineStyle',':')
% 
% ylabel('Autocorrelation','FontSize',20)
% ylim([-0.1 0.15])
% xlim([0 500])
% % set(gca,'ycolor','r')
% 
% 
% 
% 
% % decreasing neurons
% 
% % yyaxis right
% p2K0 = plot_signal_ci(lags,decACFs(:,:,1),'b-'); % Ksyn = 0
% p2K5 = plot_signal_ci(lags,decACFs(:,:,2),'b--'); % Ksyn = 0.5
% set(p2K5,'LineStyle','--')
% % p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% % set(p2K1,'LineStyle',':')
% 
% % label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% % label_h.Position(1) = 0.9;
% % label_h.Position(2) = 1;
% % set(gca,'ycolor','b')
% 
% 
% 
% % legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
% %     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
% %     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
% legend([p1K0 p1K5 p2K0 p2K5 ],...
%     'Inc Neurons, Ksyn=0.0','Inc, Ksyn=0.5',...
%     'Dec Neurons, Ksyn=0.0','Dec, Ksyn=0.5')
% 
% set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
% legend('boxoff')
% 
% grid off
% xlabel('Time Lag (ms)','FontSize',20)
% % title('Early Ramping in Slow Synapse Nets','FontSize',24)
% 
% 
% 
% % if strcmp(nettype,'flux'),ylim([-6 26]),end
% % plot_signal_ci(time/1000,EEG_ga' * 500,'g')
% saveas(gcf,['figures/' nettype '_ACFs_neurons.svg'])
% saveas(gcf,['figures/' nettype '_ACFs_neurons.png'])
% 
% 
% %% Autocorrelation Clusters
% 
% figure('Position',[0 0 450 400])
% hold on
% 
% 
% lags = 0:1000;
% 
% 
% 
% % Increasing neurons
% p1K0 = plot_signal_ci(lags,incACFs_clust(:,:,1),'r-'); % Ksyn = 0
% p1K5 = plot_signal_ci(lags,incACFs_clust(:,:,2),'r--'); % Ksyn = 0.5
% set(p1K5,'LineStyle','--')
% % p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% % set(p1K1,'LineStyle',':')
% 
% ylabel('Autocorrelation','FontSize',20)
% % ylim([-0.1 0.2])
% % xlim([0 500])
% % set(gca,'ycolor','r')
% 
% 
% 
% 
% % decreasing neurons
% 
% % yyaxis right
% p2K0 = plot_signal_ci(lags,decACFs_clust(:,:,1),'b-'); % Ksyn = 0
% p2K5 = plot_signal_ci(lags,decACFs_clust(:,:,2),'b--'); % Ksyn = 0.5
% set(p2K5,'LineStyle','--')
% % p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% % set(p2K1,'LineStyle',':')
% 
% % label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% % label_h.Position(1) = 0.9;
% % label_h.Position(2) = 1;
% % set(gca,'ycolor','b')
% 
% 
% 
% % legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
% %     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
% %     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
% legend([p1K0 p1K5 p2K0 p2K5 ],...
%     'Inc Cluster, Ksyn=0.0','Inc Cluster, Ksyn=0.5',...
%     'Other Clusters, Ksyn=0.0','Other Clusters, Ksyn=0.5')
% 
% set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
% legend('boxoff')
% 
% grid off
% xlabel('Time Lag (ms)','FontSize',20)
% % title('Early Ramping in Slow Synapse Nets','FontSize',24)
% 
% 
% 
% % if strcmp(nettype,'flux'),ylim([-6 26]),end
% % plot_signal_ci(time/1000,EEG_ga' * 500,'g')
% saveas(gcf,['figures/' nettype '_ACFs_clusters.svg'])
% saveas(gcf,['figures/' nettype '_ACFs_clusters.png'])
% 
% %%
% 
%% Autocorrelation

figure('Position',[0 0 400 320])
hold on


lags = 0:1000;

ACFsC = ACFs_clust;
ACFsN = ACFs_neuron;

ACFsC(:,lags <= 10,:) = nan;
ACFsN(:,lags <= 10,:) = nan;

lags(lags <=10) = nan;

ACFsC = ACFsC(:,~isnan(lags),:);
ACFsN = ACFsN(:,~isnan(lags),:);

lags = lags(~isnan(lags));

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

% Increasing neurons
p1K0 = plot_signal_ci(lags,ACFsN(:,:,1),'Color','k'); % Ksyn = 0
p1K5 = plot_signal_ci(lags,ACFsN(:,:,2),'Color','k'); % Ksyn = 0.5
set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
xlim([0 1000])
ylim([-0.1 0.45])
% xlim([0 500])
% set(gca,'ycolor','r')




% decreasing neurons

% yyaxis right
p2K0 = plot_signal_ci(lags,ACFsC(:,:,1),'Color','m'); % Ksyn = 0
p2K5 = plot_signal_ci(lags,ACFsC(:,:,2),'Color','m'); % Ksyn = 0.5
set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')



% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Neurons, 100% Fast','Neurons, 50% Slow',...
    'Clusters, 100% Fast','Clusters, 50% Slow')

set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_ACFs.svg'])
saveas(gcf,['figures/' nettype '_ACFs.png'])

%% Autocorrelation Neurons Only
figure('Position',[0 0 400 320])
hold on


lags = 0:1000;

ACFsC = ACFs_clust;
ACFsN = ACFs_neuron;

ACFsC(:,lags <= 10,:) = nan;
ACFsN(:,lags <= 10,:) = nan;

lags(lags <=10) = nan;

ACFsC = ACFsC(:,~isnan(lags),:);
ACFsN = ACFsN(:,~isnan(lags),:);

lags = lags(~isnan(lags));

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

% Increasing neurons
p1K0 = plot_signal_ci(lags,ACFsN(:,:,1),'Color','k'); % Ksyn = 0
p1K5 = plot_signal_ci(lags,ACFsN(:,:,2),'Color','k'); % Ksyn = 0.5
set(p1K5,'LineStyle','--')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
xlim([0 1000])
ylim([-0.1 0.45])
% xlim([0 500])
% set(gca,'ycolor','r')




% decreasing neurons

% yyaxis right
% p2K0 = plot_signal_ci(lags,ACFsC(:,:,1),'Color','m'); % Ksyn = 0
% p2K5 = plot_signal_ci(lags,ACFsC(:,:,2),'Color','m'); % Ksyn = 0.5
% set(p2K5,'LineStyle','--')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')



% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    '100% Fast','50% Slow')
set(legend,'Location','NorthEast','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_ACFsN.svg'])
saveas(gcf,['figures/' nettype '_ACFsN.png'])
%% Autocorrelation both sides

figure('Position',[0 0 400 320])
hold on

ACFsC_flip = flip(ACFs_clust,2);
ACFsN_flip = flip(ACFs_neuron,2);

ACFsC_t = cat(2,ACFsC_flip(:,1:1000,:),ACFs_clust);
ACFsN_t = cat(2,ACFsN_flip(:,1:1000,:),ACFs_neuron);

lags = -1000:1000;

ACFsC_t(:,lags >= -10 & lags <= 10,:) = nan;
ACFsN_t(:,lags >= -10 & lags <= 10,:) = nan;

lags(lags >= -10 & lags <=10) = nan;

ACFsC_t = ACFsC_t(:,~isnan(lags),:);
ACFsN_t = ACFsN_t(:,~isnan(lags),:);

lags = lags(~isnan(lags));

% ACFsC_t = sign(ACFsC_t) .* log(abs(ACFsC_t));
% ACFsN_t = sign(ACFsN_t) .* log(abs(ACFsN_t));
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

% Increasing neurons
p1K0 = plot_signal_ci(lags,ACFsN_t(:,:,1),'Color','black'); % Ksyn = 0
p1K5 = plot_signal_ci(lags,ACFsN_t(:,:,2),'Color','black'); % Ksyn = 0.5
set(p1K5,'LineStyle',':')
% p1K1 = plot_signal_ci(time/1000,incSpks_ga(:,:,3)','r:'); % Ksyn = 1
% set(p1K1,'LineStyle',':')

ylabel('Autocorrelation','FontSize',20)
% ylim([-0.1 0.45])
% xlim([0 500])
% set(gca,'ycolor','r')




% Clusters 

% yyaxis right
p2K0 = plot_signal_ci(lags,ACFsC_t(:,:,1),'Color','magenta'); % Ksyn = 0
p2K5 = plot_signal_ci(lags,ACFsC_t(:,:,2),'Color','magenta'); % Ksyn = 0.5
set(p2K5,'LineStyle',':')
% p2K1 = plot_signal_ci(time/1000,decSpks_ga(:,:,3)','b:'); % Ksyn = 1
% set(p2K1,'LineStyle',':')

% label_h = ylabel('\DeltaFiring Rate (Hz)','FontSize',20,'Rotation',270);
% label_h.Position(1) = 0.9;
% label_h.Position(2) = 1;
% set(gca,'ycolor','b')



% legend([p1K0 p1K5 p1K1 p2K0 p2K5  p2K1],...
%     'Inc, Ksyn=0.0','Inc, Ksyn=0.5','Inc, Ksyn=1',...
%     'Dec, Ksyn=0.0','Dec, Ksyn=0.5','Dec, Ksyn=1')
legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Neurons, 100% Fast','Neurons, 50% Slow',...
    'Clusters, 100% Fast','Clusters, 50% Slow')

set(legend,'Location','NorthWest','FontSize',14,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time Lag (ms)','FontSize',20)
% title('Early Ramping in Slow Synapse Nets','FontSize',24)



% if strcmp(nettype,'flux'),ylim([-6 26]),end
% plot_signal_ci(time/1000,EEG_ga' * 500,'g')
saveas(gcf,['figures/' nettype '_ACFs_flip.svg'])
saveas(gcf,['figures/' nettype '_ACFs_flip.png'])

%% Ramping Onsets

nettype = 'flux'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['./data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);

% Loading for pre-allocation
load([fpath files(1).name],'time')

time2use = (time > -3500 & time < -1000);

incOnsets = [];
decOnsets = [];

for ii=1:length(files)
    fname = files(ii).name;

    load([fpath,fname],'spks','type','time')
    
    spks_sm = squeeze(nanmean(smoothdata(spks,3,'gaussian',400) * 400,1));
    
    inc_ix = find(ismember(type,{'inc'}));
    nInc = sum(ismember(type,{'inc'}));
    
    dec_ix = find(ismember(type,{'dec'}));
    nDec = sum(ismember(type,{'dec'}));
    
    disp(['nInc: ',num2str(nInc),', nDec: ', num2str(nDec)])
    incOnsets_this = zeros(nInc,1);
    for jj = 1:nInc
    %     sigs = find(incSpks_sig(slowrampix(ii),:,1) < 0.01 & incSpks_sig(slowrampix(ii),:,2) > 0);
    %     incOnsets(ii) = time(sigs(1));
        incOnsets_this(jj) = time(findchangepts(spks_sm(inc_ix(jj),time2use),'Statistic','linear'));
    end

    decOnsets_this = zeros(nDec,1);
    for jj = 1:nDec
    %     sigs = find(incSpks_sig(slowrampix(ii),:,1) < 0.01 & incSpks_sig(slowrampix(ii),:,2) > 0);
    %     incOnsets(ii) = time(sigs(1));
        decOnsets_this(jj) = time(findchangepts(spks_sm(dec_ix(jj),time2use),'Statistic','linear'));
    end
    
    incOnsets = [incOnsets; incOnsets_this];
    decOnsets = [decOnsets; decOnsets_this];
end


figure('Position',[0 0 400 300])
histogram(incOnsets,...
    'FaceColor',[1 0.3 0.2],'LineWidth',1.5,'EdgeAlpha',1,'FaceAlpha',0.5,...
    'BinEdges',-3500:250:0,'Normalization','Probability')
hold on
histogram(decOnsets,...
    'FaceColor',[0 0.5 1],'LineWidth',1.5,'EdgeAlpha',1,'FaceAlpha',0.5,...
    'BinEdges',-3500:250:0,'Normalization','Probability')

ax = gca;
axis off
ax.YAxis.Visible = 'on';
ax.XAxis.Visible = 'on';
ax.FontSize = 16; 
ylabel('Frequency','FontSize',20) 
xlabel('Ramping Onset (s)','FontSize',20) 

xticklabels({'-3', '-2', '-1', '0'})


saveas(gcf,['figures/ramponset_sim.svg'])
saveas(gcf,['figures/ramponset_sim.png'])
    
    
