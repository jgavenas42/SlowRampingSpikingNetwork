%% Snap Exc / Inh synaptic comparisons


%% Initial
clear all
close all

addpath(genpath('./utils'))

%% Set-up

nettype = 'flux'; % flux or ignition

baseline = 500; % how many milliseconds after trial start to baseline at? 
% e.g. 500 is -3 to -2.5

fpath = ['./data/ExcInhCF/' nettype '/'];
files = dir([fpath,'*.mat']);

% Loading for pre-allocation
load([fpath files(1).name],'EEG','incSpks','decSpks',...
    'time','syntype')


incSpks_ga = nan(length(incSpks),length(files));
decSpks_ga = nan(length(decSpks),length(files));

EEG_ga = nan(length(EEG),length(files),3);

syntypes = cell(length(files),1);


for ii=1:length(files)
    fname = files(ii).name;
    load([fpath,fname],'EEG','incSpks','decSpks','syntype')
    
    
    
    EEG_ga(:,ii) = EEG - nanmean(EEG(1:baseline));
    incSpks_ga(:,ii) = incSpks - nanmean(incSpks(1:baseline));
    decSpks_ga(:,ii) = decSpks - nanmean(decSpks(1:baseline));
    
    syntypes{ii} = syntype;
    
end


incSpks_exc = incSpks_ga(:,strcmp(syntypes,'exc'));
decSpks_exc = decSpks_ga(:,strcmp(syntypes,'exc'));
EEG_exc = EEG_ga(:,strcmp(syntypes,'exc'));

incSpks_inh = incSpks_ga(:,strcmp(syntypes,'inh'));
decSpks_inh = decSpks_ga(:,strcmp(syntypes,'inh'));
EEG_inh = EEG_ga(:,strcmp(syntypes,'inh'));


%% Plotting



figure('Position',[0 0 450 400])

hold on

% Plot Excitatory slow
p1 = plot_signal_ci(time/1000,incSpks_exc(:,:)','r-');
set(p1,'Color',[0.7 0 0])

% Plot Inhibitory slow
p2 = plot_signal_ci(time/1000,decSpks_exc(:,:)','b-');
set(p2,'Color',[0 0 0.7])

xlabel('Time (s); 0 = TX','FontSize',20)
ylabel('\Delta Firing Rate (Hz)','FontSize',20)
title('Slow Excitation Only','FontSize',24)

yyaxis right
p3 = plot_signal_ci(time/1000,EEG_exc(:,:)','g-');
set(p3,'Color',[0 0.7 0])
if strcmp(nettype,'ignition')
    ylim([-0.02 0.06])
else
    ylim([-0.02 0.05])
end

set(gca,'ycolor',[0 0.7 0]) 
label_h = ylabel('Voltage (\muV)','FontSize',20,'Rotation',270);

label_h.Position(1) = 1;
% label_h.Position(2) = 1;





p4 = yline(0,'k--');
p5 = xline(0,'k--');

legend([p1 p2 p3],...
    'Increasing Neurons','Decreasing Neurons', 'EEG')

set(legend,'Location','NorthWest','FontSize',16)

legend('boxoff')

grid off

saveas(gcf,['figures/' nettype '_ExcInhCF_exc.svg'])

%% Slow Inhibition

figure('Position',[0 0 450 400])
hold on

% Plot Excitatory slow
p1 = plot_signal_ci(time/1000,incSpks_inh(:,:)','r-');
set(p1,'Color',[0.7 0 0])

% Plot Inhibitory slow
p2 = plot_signal_ci(time/1000,decSpks_inh(:,:)','b-');
set(p2,'Color',[0 0 0.7])

xlabel('Time (s); 0 = TX','FontSize',20)
ylabel('\DeltaFiring Rate (Hz)','FontSize',20)
title('Slow Inhibition Only','FontSize',24)

yyaxis right
p3 = plot_signal_ci(time/1000,EEG_inh(:,:)','g-');
set(p3,'Color',[0 0.7 0])

if strcmp(nettype,'ignition')
    ylim([-0.15 0.3])
else
    ylim([-0.14 0.7])
end
set(gca,'ycolor',[0 0.7 0]) 

label_h = ylabel('Voltage (\muV)','FontSize',20,'Rotation',270);

label_h.Position(1) = 1;



% label_h.Position(2) = 0.25;



p4 = yline(0,'k--');
p5 = xline(0,'k--');

legend([p1 p2 p3],...
    'Increasing Neurons','Decreasing Neurons', 'EEG')

set(legend,'Location','NorthWest','FontSize',16)

legend('boxoff')

grid off

saveas(gcf,['figures/' nettype '_ExcInhCF_inh.svg'])

