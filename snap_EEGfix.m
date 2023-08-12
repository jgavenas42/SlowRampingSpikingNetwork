%% SNAP EEG Fix

% Compare EEG proxy with all PSPs vs only PSPs onto excitatory neurons

nettype = 'flux_eegfix'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['./data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);

% Loading for pre-allocation
load([fpath files(1).name],'EEG','EEG_inh','time')


EEG_ga = nan(length(EEG),length(files),3);
EEG_fix_ga = nan(length(EEG_inh),length(files),3);




fpath = ['./data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);


nInc = nan(length(files),3);
nDec = nan(length(files),3);
nNon = nan(length(files),3);


for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','EEG_inh','time','tx')
    
    if sum(~isnan(tx)) < 15, continue, end
    
    
    EEG_ga(:,ii,1) = EEG - nanmean(EEG(1:baseline));
    EEG_fix_ga(:,ii,1) = EEG_inh - nanmean(EEG_inh(1:baseline));
    
    
    
end


% Add in the Ksyn 0.5
fpath = ['./data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);

for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','EEG_inh','time','tx')
    
    if sum(~isnan(tx)) < 15, continue, end
    
    
    EEG_ga(:,ii,2) = EEG - nanmean(EEG(1:baseline));
    EEG_fix_ga(:,ii,2) = EEG_inh - nanmean(EEG_inh(1:baseline));
    
    
    
end

%% EEG Regular
figure('Position',[0 0 800 320])

tiledlayout(1,2)

ax1 = nexttile;

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

title('PSPs onto exc neurons only','FontSize',24)


grid off
xlabel('Time rel. threshold-crossing (s)','FontSize',20)
ylabel('Voltage (mV)','FontSize',20)

%% EEG fixed -- PSPs onto all neurons
% ax2=nexttile;

figure('Position',[0 0 400 320])

hold on
p1 = plot_signal_ci(time/1000,EEG_fix_ga(:,:,1)','Color',[0 0.5 0]);
p2 = plot_signal_ci(time/1000,EEG_fix_ga(:,:,2)','Color',[0 0.5 0]);
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

% title('PSPs onto all neurons','FontSize',24)

% linkaxes([ax1 ax2],'y')

ylim([-0.02,0.005])

saveas(gcf,['figures/' nettype '_EEG.svg'])
saveas(gcf,['figures/' nettype '_EEG.png'])
