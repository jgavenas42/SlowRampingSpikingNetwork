%% Quantifying the Charactieristics of Spontaneous Fluctuations

% By quantifying the presence of fluctuations using CV, and the presence
% of bistability using the ST measure (Schaub et al., 2015), we will
% compare the effects of clustering on fluctuations in the presence of fast
% and slow synapses

%% Testing set-up
Ree_list = 1:0.1:6;
connMat_list = 1:50;
Ksyn_list = 0:0.1:0.5;

N = 400;
Nclusters = 4;
noisescale = 0.5;
Jee = 0.2;
Jclust = 1.9;
tsteps = 5000;

table_length = length(Ree_list) *length(connMat_list) * length(Ksyn_list);
tab_Ree = nan(table_length,1);
tab_rng = nan(table_length,1);
tab_ksyn = nan(table_length,1);
tab_cv = nan(table_length,1);
tab_st = nan(table_length,1);
tab_ss = nan(table_length,1);
tab_ff = nan(table_length,1);
tab_fexp = nan(table_length,1);
tab_tau = nan(table_length,1);

powspectrum = nan(table_length,89); % have to dig in a bit
acfs = nan(table_length,20);

counter = 1;
tic
for ii =1:length(Ree_list)
    for jj = 1:length(connMat_list)
        for kk = 1:length(Ksyn_list)
            
            Ree = Ree_list(ii);
            connMat_rng = connMat_list(jj);
            Ksyn = Ksyn_list(kk);




            % Create the connectivity matrix
            W = connMatrix('N',N,...
                'Nclusters',Nclusters,...
                'Ree',Ree,...
                'Jee',Jee,...
                'Jclust',Jclust,...
                'connMat_rng',connMat_rng);


            [spks,Rs] = snap_singletrial(W,'tsteps',tsteps,...
                            'noisescale',noisescale,...
                            'Ksyn',Ksyn); 

            ST = calculateST(spks(1:320,501:end),Nclusters);
            CV = calculateCV(spks(1:320,501:end));
%             SS = calculateSilhouetteScore(spks(1:320,501:end),2);
%             FF = calculateFanoFactor(spks(1:320,501:end));
            [fexp,Pxx,frex] = calculatePowerSpectrum(spks(1:320,501:end),0);
            [tau,acf,lags] = calculateACFtau(spks(1:320,501:end),1);

            saveas(gcf,sprintf('../figures/flux_powspectrums/%0.2f_%0.2f_%d.png',Ree,Ksyn,connMat_rng))
            close

            tab_Ree(counter) = Ree;
            tab_rng(counter) = connMat_rng;
            tab_ksyn(counter) = Ksyn;
            tab_cv(counter) = CV;
            tab_st(counter) = ST;
            tab_fexp(counter) = fexp;
            tab_tau(counter) = tau;
            
            powspectrum(counter,:) = Pxx;
            acfs(counter,:) = acf;

            counter = counter + 1;
            fprintf('Ree = %0.2f, Network = %d, Ksyn = %0.2f, 1/f exp = %0.2f, elapsed time = %0.2f\n',Ree,connMat_rng,Ksyn,fexp,toc)
            
        end
    end
end


fluxQuant = table;

fluxQuant.Ree = tab_Ree;
fluxQuant.Network = tab_rng;
fluxQuant.Ksyn = tab_ksyn;
fluxQuant.CV = tab_cv;
fluxQuant.ST = tab_st;
fluxQuant.fexp = tab_fexp;
fluxQuant.taus = tab_tau;

writetable(fluxQuant,'../fluxQuant/fluxquant.csv')

data = [];
data.Ree = tab_Ree;
data.Network = tab_rng;
data.Ksyn = tab_ksyn;
data.frex = frex;
data.powspectrum = powspectrum;
data.lags = lags;
data.acf = acfs;
save('../fluxQuant/powspectrums.mat','data')

%% Plotting

% Plot power spectrums for 

Ree_quintiles = [1; 2; 3; 4; 5; 6];

Ree_ix = [];

for ii = 1:length(Ree_quintiles)-1
    Ree_ix = [Ree_ix (data.Ree >= Ree_quintiles(ii) & data.Ree < Ree_quintiles(ii+1))];
end

Ksyn_ix = [];

for ii = 1:length(Ksyn_list)
    Ksyn_ix = [Ksyn_ix data.Ksyn == Ksyn_list(ii)];
end

powspectrums_avg = [];

Ree2use = [];
Ksyn2use = [];
for ii = 1:length(Ksyn_list)
    for jj = 1:length(Ree_quintiles)-1
        powspectrums_avg = [powspectrums_avg; nanmean(data.powspectrum(Ree_ix(:,jj) == 1 & Ksyn_ix(:,ii) == 1,:),1)];
        Ree2use = [Ree2use Ree_quintiles(jj)];
        Ksyn2use = [Ksyn2use Ksyn_list(ii)];
    end
end

figure('Position',[0 0 1400 320])
hold on

cmap = brewermap(ceil(length(Ksyn_list)*1.2),'-RdYlBu');

cmap = [0.4176 0.0006 0.6584;...
         0.5897 0.0729 0.6304;...
         0.7360 0.2094 0.5279;...
         0.8501 0.3470 0.4172;...
         0.9380 0.4927 0.3126;...
         0.9883 0.6523 0.2114];

tileplot = tiledlayout(1,length(Ree_quintiles));
tileplot = tiledlayout('flow');


for jj = 1:length(Ree_quintiles)-1
    
    nexttile
    hold on
    for ii = 1:length(Ksyn_list)
        plot(data.frex,powspectrums_avg(Ree2use == Ree_quintiles(jj) & Ksyn2use == Ksyn_list(ii),:),...
            'Color',cmap(ii,:),'LineWidth',2.5)
%         plot(data.frex,log(powspectrums_avg(Ree2use == Ree_quintiles(jj) & Ksyn2use == Ksyn_list(ii),:)),...
%             'Color',cmap(ii+1,:),'LineWidth',2.5)
    end
    
    set(gca, 'YScale', 'log','XScale','log')
    
    alpha 0.7
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('Frequency (Hz)','FontSize',14)
    xticks([1 5 10])
    
    if jj < 5
        title(sprintf("R_{EE} \\in [%0.2f,%0.2f)",Ree_quintiles(jj),Ree_quintiles(jj+1)),'FontSize',16,'Interpreter','tex')
    else
        title(sprintf("R_{EE} \\in [%0.2f,%0.2f]",Ree_quintiles(jj),Ree_quintiles(jj+1)),'FontSize',16,'Interpreter','tex')
    end
    
    
%     ylim([min(log(powspectrums_avg),[],'all')-0.5 max(log(powspectrums_avg),[],'all')+0.5])
    ylim([min(powspectrums_avg,[],'all')-10 max(powspectrums_avg,[],'all')+1000])
    
    if jj == 1
        ylabel('Spike-Train Power (1/s)','FontSize',14)

    else
        ax.YTick = [];
    end
    

    
end

% colormap(cmap(2:7,:));
colormap(cmap(1:6,:));

a = colorbar;
a.Ticks = linspace(0.04,0.46,6);
a.TickLabels = {'0%','10%','20%','30%','40%','50%'};
a.TickLength = 0.;

a.FontSize = 14;
caxis([0 0.5])
ylabel(a,'% Slow Synapses','FontSize',16,'Rotation',270);
a.Label.Position(1) = 8;
a.Label.Position(2) = 0.25;

tileplot.TileSpacing = 'compact';
% tileplot.Padding = 'compact';

saveas(gcf,['../figures/fluctuations/powspectrum_loglog.svg'])

saveas(gcf,['../figures/fluctuations/powspectrum_loglog.png'])

%% Autocorrelation Plotting
acf_avg = [];

Ree2use = [];
Ksyn2use = [];
for ii = 1:length(Ksyn_list)
    for jj = 1:length(Ree_quintiles)-1
        acf_avg = [acf_avg; nanmean(data.acf(Ree_ix(:,jj) == 1 & Ksyn_ix(:,ii) == 1,:),1)];
        Ree2use = [Ree2use Ree_quintiles(jj)];
        Ksyn2use = [Ksyn2use Ksyn_list(ii)];
    end
end

figure('Position',[0 0 1400 320])
hold on

% cmap = brewermap(ceil(length(Ksyn_list)*1.2),'-RdYlBu');
cmap = [0.4176 0.0006 0.6584;...
         0.5897 0.0729 0.6304;...
         0.7360 0.2094 0.5279;...
         0.8501 0.3470 0.4172;...
         0.9380 0.4927 0.3126;...
         0.9883 0.6523 0.2114];
tileplot = tiledlayout(1,length(Ree_quintiles));
tileplot = tiledlayout('flow');


for jj = 1:length(Ree_quintiles)-1
    
    nexttile
    hold on
    for ii = 1:length(Ksyn_list)
        plot(data.lags,acf_avg(Ree2use == Ree_quintiles(jj) & Ksyn2use == Ksyn_list(ii),:),...
            'Color',cmap(ii,:),'LineWidth',2.5)
%         plot(data.frex,log(powspectrums_avg(Ree2use == Ree_quintiles(jj) & Ksyn2use == Ksyn_list(ii),:)),...
%             'Color',cmap(ii+1,:),'LineWidth',2.5)
    end
    
%     set(gca, 'YScale', 'log')
    
    alpha 0.7
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('Time Lag (ms)','FontSize',14)
    
    if jj < 5
        title(sprintf("R_{EE} \\in [%0.2f,%0.2f)",Ree_quintiles(jj),Ree_quintiles(jj+1)),'FontSize',16,'Interpreter','tex')
    else
        title(sprintf("R_{EE} \\in [%0.2f,%0.2f]",Ree_quintiles(jj),Ree_quintiles(jj+1)),'FontSize',16,'Interpreter','tex')
    end
    
    
%     ylim([min(log(powspectrums_avg),[],'all')-0.5 max(log(powspectrums_avg),[],'all')+0.5])
    ylim([min(acf_avg,[],'all')-0.01 max(acf_avg,[],'all')+0.01])
    
    if jj == 1
        ylabel('Autocorrelation','FontSize',14)

    else
        ax.YTick = [];
    end
    

    
end

colormap(cmap(1:6,:));
% colormap(cmap(2:7,:));

a = colorbar;
a.Ticks = linspace(0.04,0.46,6);
a.TickLabels = {'0%','10%','20%','30%','40%','50%'};
a.TickLength = 0.;

a.FontSize = 14;
caxis([0 0.5])
ylabel(a,'% Slow Synapses','FontSize',16,'Rotation',270);
a.Label.Position(1) = 8;
a.Label.Position(2) = 0.25;

tileplot.TileSpacing = 'compact';
% tileplot.Padding = 'compact';

saveas(gcf,['../figures/fluctuations/autocorrelation.svg'])
saveas(gcf,['../figures/fluctuations/autocorrelation.png'])
