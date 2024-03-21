%% SNAP goodness-of-fit

% Here we will first load in the FRs over time for real data, subselect the
% neurons that exhibit slow ramping, then loop over simulated networks and
% measure the goodness-of-fit of the networks with Ksyn 0.5 vs 0
clear all

%% Load the real firing rates
load('../data/realFRs.mat')

incSpks_real = incSpks;
decSpks_real = decSpks;
time_real = time;


%% Initial plotting


figure('Position',[0 0 600 480])
hold on
yline(0,'-k','LineWidth',1)

ymax = max(incSpks_real,[],'all');

gain_cutoff = 0.4;

p1K0 = plot_signal_ci(time_real,incSpks_real(incGain <= gain_cutoff & ismember(incRegion,mfcRegions),:),'Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(time_real,incSpks_real(incGain >= gain_cutoff & ismember(incRegion,mfcRegions),:),'Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K0,'LineStyle',':')
set(p1K5,'LineStyle','-')


yyaxis left

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis(1).FontSize = 14;
ax.YAxis(2).FontSize = 14;

ylabel('\DeltaFiring Rate (Hz)','FontSize',20)


% Setting ylimit so they work together
yl = ylim;
ylim([yl(1) - 2, yl(2) + 1])
yl = ylim;

set(gca,'ycolor',[1 0.3 0.2])

% decreasing neurons

yyaxis right

p2K0 = plot_signal_ci(time_real,decSpks_real(decGain <= gain_cutoff & ismember(decRegion,mfcRegions),:),'Color',[0 0.5 1]); % Ksyn = 0
p2K5 = plot_signal_ci(time_real,decSpks_real(decGain >= gain_cutoff & ismember(decRegion,mfcRegions),:),'Color',[0 0.5 1]); % Ksyn = 0.5
set(p2K0,'LineStyle',':')
set(p2K5,'LineStyle','-')


ylim(yl/2)

% ax.YColor = [0 0.5 1];
set(gca,'ycolor',[0 0.5 1])


legend([p1K0 p1K5 p2K0 p2K5 ],...
    'Inc, Steplike','Inc, Slow Ramping',...
    'Dec, Steplike','Dec, Slow Ramping')

set(legend,'Location','NorthWest','FontSize',16,'NumColumns',1)
legend('boxoff')

grid off
xlabel('Time rel. movement (s)','FontSize',20)
xlim([-3 0])



saveas(gcf,['../figures/real_FRs_both_',num2str(gain_cutoff),'_mfc_new.svg'])
saveas(gcf,['../figures/real_FRs_both_',num2str(gain_cutoff),'_mfc_new.png'])


%% Set-up

nettype = 'flux'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['../data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);
files(strncmp({files.name}, '.', 1)) = []; %remove files and dir starting with '.'

% Loading for pre-allocation
load([fpath files(2).name],'EEG','incFRs','decFRs','tx',...
    'incSpks','decSpks','nonSpks','time','time_fr','taus')

decFRs_ga = nan(length(decFRs),length(files),3);
incFRs_ga = nan(length(incFRs),length(files),3);

incSpks_ga = nan(length(incSpks),length(files),3);
decSpks_ga = nan(length(decSpks),length(files),3);
nonSpks_ga = nan(length(nonSpks),length(files),3);

EEG_ga = nan(length(EEG),length(files),3);

[f,xi] = ksdensity(tx/1000,0:0.5:10);

tx_ga = nan(length(f),length(files),3);


fpath = ['../data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);
files(strncmp({files.name}, '.', 1)) = []; %remove files and dir starting with '.'


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
fpath = ['../data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);
files(strncmp({files.name}, '.', 1)) = []; %remove files and dir starting with '.'


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
fpath = ['../data/' nettype '/Ksyn_' num2str(1) '/'];
files = dir([fpath,'*.mat']);
files(strncmp({files.name}, '.', 1)) = []; %remove files and dir starting with '.'

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

%% Plot

figure('Position',[0 0 600 480])
hold on
yline(0,'-k','LineWidth',1)

ymax = max(incSpks_ga,[],'all');

% Increasing neurons
p1K0 = plot_signal_ci(time/1000,incSpks_ga(:,:,1)','Color',[1 0.3 0.2]); % Ksyn = 0
p1K5 = plot_signal_ci(time/1000,incSpks_ga(:,:,2)','Color',[1 0.3 0.2]); % Ksyn = 0.5
set(p1K0,'LineStyle',':')
set(p1K5,'LineStyle','-')

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
set(p2K0,'LineStyle',':')
set(p2K5,'LineStyle','-')

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

saveas(gcf,['../figures/sim_FRs_both_new.png'])


%%

YYinc = nanmean(incSpks_real(incGain <= gain_cutoff & ismember(incRegion,mfcRegions),time_real > -3),1);
YYdec = nanmean(decSpks_real(decGain <= gain_cutoff & ismember(decRegion,mfcRegions),time_real > -3),1);

Rs = zeros(size(incSpks_ga,2),4);

for ii = 1:size(incSpks_ga,2)
    XXinc_0 = incSpks_ga(time < 0,ii,1);
    XXdec_0 = decSpks_ga(time < 0,ii,1);
    XXinc_5 = incSpks_ga(time < 0,ii,2);
    XXdec_5 = decSpks_ga(time < 0,ii,2);
    
    % increasing Ksyn 0
    mdl = fitlm(XXinc_0,YYinc);
    Rs(ii,1) = mdl.Rsquared.Ordinary;
    % decreasing Ksyn 0
    mdl = fitlm(XXdec_0,YYdec);
    Rs(ii,2) = mdl.Rsquared.Ordinary;
    % increasing Ksyn 5
    mdl = fitlm(XXinc_5,YYinc);
    Rs(ii,3) = mdl.Rsquared.Ordinary;
    % decreasing Ksyn 5
    mdl = fitlm(XXdec_5,YYdec);
    Rs(ii,4) = mdl.Rsquared.Ordinary;
end

Rs_mean = nanmean(Rs,1);
Rs_std = nanstd(Rs,1);

dof = size(Rs,1) - sum(isnan(Rs),1);
t_low = tinv([0.025 0.025 0.025 0.025],dof);
t_high = tinv([0.975 0.975 0.975 0.975],dof);

CIdist = t_high .* Rs_std ./ sqrt(size(Rs,1));


%% Plotting GOF

figure('Position',[100 100 300 580])

hold on
ax = gca;
ax.FontSize = 16; 
set(ax,'YAxisLocation','right')
colors2use = [1 0.3 0.2;0 0.5 1;1 0.3 0.2;0 0.5 1];

xpositions = [1 1.6 3.4 4];

for ii = 1:4
%     plot(randn(size(Rs,1),1)*0.1+ii,Rs(:,ii),'.','Color',colors2use(ii,:),'MarkerSize',12);
%     scatter(randn(size(Rs,1),1)*0.1+xpositions(ii),Rs(:,ii),25,colors2use(ii,:),"filled",'MarkerFaceAlpha',0.5);
    scatter((1:size(Rs,1))*0.01+xpositions(ii)-0.1,Rs(:,ii),30,colors2use(ii,:),"filled",'MarkerFaceAlpha',0.5);

    errorbar(xpositions(ii),Rs_mean(ii),CIdist(ii),'vert','.','Color','black','LineWidth',2);
end

set(gca,'TickLength',[0, 0])
xlim([0.5 4.5])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
label_h = ylabel('R^{2}','FontSize',20,'rotation',0,'HorizontalAlignment','left');
% label_h.Position(1) = 10; % change horizontal position of ylabel
% label_h.Position(2) = 0.5; % change vertical position of ylabel
% label_h = ylabel('R^{2} Simulated (50% slow) vs Real (Slow Ramping)','FontSize',20,'rotation',90,'HorizontalAlignment','left');

xticklabels(["100% Fast","50% slow"]);
xtickangle(ax,30)

% set(ax,'Position', [20 0 180 480]);

saveas(gcf,['../figures/sims_rsquared.png'])

%%
% stats tests

signrank(Rs(:,1),Rs(:,3))



