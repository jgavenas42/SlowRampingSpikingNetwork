%% Shared neural correlations
clear all 
close all

addpath(genpath('./utils'))

% Want to test pair-wise correlations between neurons and 



%% Testing set-up
Ree = 3.4;
N = 400;
Nclusters = 4;
noisescale = 0.5;
Ksyn = 0.5;
% Testing
% Ksyn = [randn(320,1)*0.1+0.5; ones(80,1)*0.];
% Ksyn = [ones(320,1)*0.; ones(80,1)*0.5];
% Ksyn(1:2:end) = 0;

connMat_rng = 5;
Jee = 0.2;
Jclust = 1.9;

tsteps = 5000;


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
            
%% Plot Example     
% figure 
figure('Position',[0 0 550 500])
tiledlayout(3,1,'TileSpacing','compact')
nexttile([2 1])
plotRaster(spks,'MarkerSize',5,'ColorSep',true,'Color_exc',[0.25 0 0.25],'Color_inh',[0 0.5 0]);
ylim([0,size(spks(~isnan(spks(:,1,1)),:,:),1)+1])
set(gca,'FontSize',24,'XTick',[])
% xlabel('Time (s)','FontSize',24)
% ylabel('Neuron Index','FontSize',24)
% xticks(0:1000:5000)
% xticklabels(0:1:5)

% saveas(gcf,['figures/examples/rasterEx_Ree_',num2str(Ree),'_Ksyn_',num2str(Ksyn), '_Net_',num2str(connMat_rng),'.svg'])
% saveas(gcf,['figures/examples/rasterEx_Ree_',num2str(Ree),'_Ksyn_',num2str(Ksyn), '_Net_',num2str(connMat_rng),'.png'])

% figure
% FRs = firingRate(spks);
% cluster_FRs = squeeze(mean(reshape(FRs(1:0.8*N,:),0.8*N/Nclusters,Nclusters,[]),1));
% plot((1:50:tsteps)/(1000),cluster_FRs,'LineWidth',2)
% set(gca,'FontSize',16)
% xlabel('Time (s)','FontSize',20)
% ylabel('Firing Rate (Hz)','FontSize',20)
% legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'},...
%     'Location','northwest')
% legend boxoff

nexttile
spks_sm = smoothdata(spks,2,'gaussian',400);
cluster_spks = squeeze(mean(reshape(spks_sm(1:0.8*N,:),0.8*N/Nclusters,Nclusters,[]),1));
plot((1:tsteps)/(1000),cluster_spks*400,'LineWidth',2)
set(gca,'FontSize',24)
xlabel('Time (s)','FontSize',32)
ylabel('Firing Rate (Hz)','FontSize',28)
% legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'},...
%     'Location','northwest')
% legend boxoff

saveas(gcf,['./figures/examples/Net_',num2str(connMat_rng),'_Ksyn_',num2str(Ksyn), '_Ree_',num2str(Ree),'.svg'])
saveas(gcf,['./figures/examples/Net_',num2str(connMat_rng),'_Ksyn_',num2str(Ksyn), '_Ree_',num2str(Ree),'.png'])



% disp(['ST = ',num2str(calculateST(spks(1:320,501:end),4)),...
%     ', CV = ',num2str(calculateCV(spks(1:320,501:end))),...
%     ', FF = ',num2str(calculateFanoFactor(spks(1:320,501:end))),...
%     ', Silhouette Score = ',num2str(calculateSilhouetteScore(spks(1:40:320,501:end),2))])
%% Testing
cluster_spks = squeeze(mean(reshape(spks_sm(:,:),0.8*N/Nclusters,Nclusters+1,[]),1));

cluster_spks_zoom = zeros(3,size(cluster_spks,2));
[t_, max_ix] = max(mean(cluster_spks(1:Nclusters,:),2));

cluster_spks_zoom(1,:) = cluster_spks(max_ix,:);
cluster_spks_zoom(2,:) = mean(cluster_spks(setdiff(1:end,[max_ix 5]),:),1);
cluster_spks_zoom(3,:) = cluster_spks(end,:);

figure
plot((1:tsteps)/(1000),cluster_spks_zoom(1,:)*400,...
    'Color','red','LineWidth',2);
hold on
plot((1:tsteps)/(1000),cluster_spks_zoom(2,:)*400,...
    'Color','#EDB120','LineWidth',2);
plot((1:tsteps)/(1000),cluster_spks_zoom(3,:)*400,...
    'Color','b','LineWidth',2);

set(gca,'FontSize',16)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',20)
legend('Active Cluster','Other Clusters (Avg)','Inhibitory Pop.',...
    'Location','northwest')
legend boxoff
% xlim([6 10])
    
% saveas(gcf,['./final_figures/clusterEx_zoom_Ree_',num2str(Ree),'_Ksyn_',num2str(Ksyn), '_Net_',num2str(connMat_rng),'.epsc'])
% saveas(gcf,['./final_figures/clusterEx_zoom_Ree_',num2str(Ree),'_Ksyn_',num2str(Ksyn), '_Net_',num2str(connMat_rng),'.png'])
%%
figure
FRs = firingRate(spks);
spks_sm = smoothdata(spks,2,'gaussian',400);
cluster_spks = squeeze(mean(reshape(spks_sm(1:0.8*N,:),0.8*N/Nclusters,Nclusters,[]),1));
plot((1:tsteps)/(1000),cluster_spks,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',20)
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'},...
    'Location','northwest')
legend boxoff

%% Calculating Correlations
nCk = nchoosek(1:N*0.8,2);

noise_corr = zeros(length(nCk),1);
shareclust = zeros(length(nCk),1);
spks_sm = smoothdata(spks,2,'gaussian',400);
% spks_sm = firingRate(spks);

tic
for ii = 1:length(nCk)
    
    % Get 2x2 matrix of correlation coefficients
    corr_temp = corr([spks_sm(nCk(ii,1),:)' spks_sm(nCk(ii,2),:)'],...
        'Type','Spearman');
    
    noise_corr(ii) = corr_temp(1,2); % Get correlation of interest
    
    if nCk(ii,1) < 81 && nCk(ii,2) < 81 
        shareclust(ii)=1;
    elseif (nCk(ii,1)>=81 && nCk(ii,1) < 161) && (nCk(ii,2)>=81 && nCk(ii,2) < 161)
        shareclust(ii)=2;
    elseif (nCk(ii,1)>=161 && nCk(ii,1) < 241) && (nCk(ii,2)>=161 && nCk(ii,2) < 241)
        shareclust(ii)=3;
    elseif (nCk(ii,1)>=241 && nCk(ii,1) < 321) && (nCk(ii,2)>=241 && nCk(ii,2) < 321)
        shareclust(ii)=4;
    end
    if mod(ii,1000) == 1,toc,end
end


%%

nc_stats = zeros(5,3);

for ii=0:Nclusters
    x = noise_corr((~isnan(noise_corr))&(shareclust==ii));
    
    SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    CI = mean(x) + ts*SEM;
    
    nc_stats(ii+1,1) = mean(x);
    nc_stats(ii+1,2) = CI(1);
    nc_stats(ii+1,3) = CI(2);
end
    
X = categorical({'Unshared','C1','C2','C3','C4'});
X = reordercats(X,{'Unshared','C1','C2','C3','C4'});

Y = nc_stats(:,1);
% Make a bar plot
figure
bar(X,Y)

hold on

errorbar([1 2 3 4 5],[nc_stats(:,1)],[nc_stats(:,3)-nc_stats(:,1)],'o')



%% Real data

fpath = ('data_real/');

filenames = dir([fpath,'*.mat']);


nois = 1:36; % sub 10 ses 1
% nois = 37:54; % sub 10 ses 2 
% nois = 55:111; % sub 11 ses 1
% nois = 112:121; % sub 12 ses 2 no results (<2 inc and dec)
% nois = 122:159; % sub 1 ses 1 less than dec-dec
% nois = 160:217; % sub 2 ses 1 only one inc-inc pair
% nois = 218:273; % sub 2 ses 2
% nois = 274:314; % sub 3 ses 1 no results (<2 inc)
% nois = 315:352; % sub 4 ses 1 no results (<2 inc)
% nois = 353:394; % Sub 4 ses 2
% nois = 395:433; % sub 4 ses 3 negative correlation
% nois = 434:482; % sub 4 ses 4 only 2 neurons
% nois = 483:525; % sub 4 ses 5
% nois = 526:569; % sub 4 ses 6 no results (<2 inc)
% nois = 570:590; % sub 5 ses 1 only 2 neurons, less than dec
% nois = 591:614; % sub 5 ses 2 no results any (<2 inc and dec)
% nois = 615:627; % sub 5 ses 3 no results (<2 inc)
% nois = 628:649; % sub 5 ses 4 only 2 neurons
% nois = 650:697; % sub 6 ses 1 no results (<2 inc)
% nois = 698:746; % sub 6 ses 2 smaller than dec
% nois = 747:773; % sub 6 ses 3 no results (<2 inc)
% nois = 774:831; % sub 6 ses 4 smaller than dec
% nois = 832:868; % sub 7 ses 1 no results
% nois = 869:895; % sub 7 ses 2 negative correlation
% nois = 896:938; % sub 8 ses 1 no results (<2 inc)
% nois = 939:995; %sub 8 ses 2 
% nois = 996:1012; % sub 9 ses 1 no results (<2 inc)
% nois = 1013:1032; % sub 9 ses 2 no results (<2 inc)


load([fpath,filenames(nois(1)).name]);

figures_path = 'figures/scratch';

nNeurons = length(nois);
spks_fmk = nan(nNeurons,1001*ntrials);
types = cell(nNeurons,1);
regions = cell(nNeurons,1);
channels = cell(nNeurons,1);
clusters = cell(nNeurons,1); 
gains = zeros(nNeurons,1);

% for ii=1:length(filenames)
for ii = 1:length(nois)
    load([fpath,filenames(ii+nois(1)-1).name],'type','BaselineAlignedSpikes',...
        'region');
    if strcmp(type,'nochange') || strcmp(type,'bad')
        gains(ii) = nan;
    else
        load([fpath,filenames(ii+nois(1)-1).name],'gain')
        gains(ii) = gain;
    end
    smdat = smoothdata(BaselineAlignedSpikes,2,'gaussian',400);
    spks_fmk(ii,:) = reshape(smdat(1:ntrials,:)',[],1);
    types{ii} = type;
    regions{ii} = region;
    channels{ii} = channel;
    clusters{ii} = cluster;
end
    
nCk = nchoosek(1:length(types),2);

spks_sm = smoothdata(spks_fmk,2,'gaussian',400);

nc_fmk = nan(length(nCk),1);
typematch = nan(length(nCk),1);
type1 = cell(length(nCk),1);
type2 = cell(length(nCk),1);
reg1 = cell(length(nCk),1);
reg2 = cell(length(nCk),1);
gain1 = nan(length(nCk),1);
gain2 = nan(length(nCk),1);


tic
for ii = 1:length(nCk)
    
    % Get 2x2 matrix of correlation coefficients
    corr_temp = corr([spks_sm(nCk(ii,1),:)' spks_sm(nCk(ii,2),:)'],...
        'Type','Spearman');
    
    nc_fmk(ii) = corr_temp(1,2); % Get correlation of interest
    
    if strcmp(types{nCk(ii,1)},'inc') & strcmp(types{nCk(ii,2)},'inc')
        typematch(ii) = 1;
    elseif strcmp(types{nCk(ii,1)},'dec') & strcmp(types{nCk(ii,2)},'dec')
        typematch(ii) = 2;
    end
    if mod(ii,1000)==0,toc,end   
    
    type1{ii} = types{nCk(ii,1)};
    type2{ii} = types{nCk(ii,2)};
    
    reg1{ii} = regions{nCk(ii,1)};
    reg2{ii} = regions{nCk(ii,2)};
    
    gain1{ii} = gains(nCk(ii,1));
    gain2{ii} = gains(nCk(ii,1));
end


nc_stats = zeros(3,3);

for ii=0:2
    if ii == 0
        x = nc_fmk(~isnan(nc_fmk));
    else
        x = nc_fmk(~isnan(nc_fmk)&typematch==ii);
    end
    
    SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    CI = mean(x) + ts*SEM;
    
    nc_stats(ii+1,1) = mean(x);
    nc_stats(ii+1,2) = CI(1);
    nc_stats(ii+1,3) = CI(2);
end


X = categorical({'Inc-Dec','Inc-Inc','Dec-Dec'});
X = reordercats(X,{'Inc-Dec','Inc-Inc','Dec-Dec'});

Y = nc_stats(:,1);
% Make a bar plot
figure
bar(X,Y)

hold on

errorbar([1 2 3],[nc_stats(:,1)],[nc_stats(:,3)-nc_stats(:,1)],'o')

pwCorrs = table(nc_fmk,typematch,type1,type2,reg1,reg2);
pwCorrs.subject = ones(size(nc_fmk))*subj;
pwCorrs.session = ones(size(nc_fmk))*session;

writetable(pwCorrs,['./pwCorr/sub',num2str(subj),'ses',num2str(session),'.csv'])



%% Done in a loop so we can try out different measures of correlation

fpath = ('data_real/');

filenames = dir([fpath,'*.mat']);


% datalist_ix gives the starting and stopping indexes for each subjects'
% data
data_ixs = [1 36; % sub 10 ses 1
    37 54; % sub 10 ses 2 
    55 111; % sub 11 ses 1
    112 121; % sub 12 ses 2 no results (<2 inc and dec)
    122 159; % sub 1 ses 1 less than dec-dec
    160 217; % sub 2 ses 1 only one inc-inc pair
    218 273; % sub 2 ses 2
    274 314; % sub 3 ses 1 no results (<2 inc)
    315 352; % sub 4 ses 1 no results (<2 inc)
    353 394; % Sub 4 ses 2
    395 433; % sub 4 ses 3 negative correlation
    434 482; % sub 4 ses 4 only 2 neurons
    483 525; % sub 4 ses 5
    526 569; % sub 4 ses 6 no results (<2 inc)
    570 590; % sub 5 ses 1 only 2 neurons, less than dec
    591 614; % sub 5 ses 2 no results any (<2 inc and dec)
    615 627; % sub 5 ses 3 no results (<2 inc)
    628 649; % sub 5 ses 4 only 2 neurons
    650 697; % sub 6 ses 1 no results (<2 inc)
    698 746; % sub 6 ses 2 smaller than dec
    747 773; % sub 6 ses 3 no results (<2 inc)
    774 831; % sub 6 ses 4 smaller than dec
    832 868; % sub 7 ses 1 no results
    869 895; % sub 7 ses 2 negative correlation
    896 938; % sub 8 ses 1 no results (<2 inc)
    939 995; %sub 8 ses 2 
    996 1012; % sub 9 ses 1 no results (<2 inc)
    1013 1032]; % sub 9 ses 2 no results (<2 inc)


for jj = 1:length(data_ixs)
    
    % Create list of data corresponding to given subject & session
    nois = data_ixs(jj,1):data_ixs(jj,2);
    
    % Load data from first neuron in that session to pre-allocate info
    load([fpath,filenames(nois(1)).name]);

    % Figures path
    figures_path = 'figures/scratch';

    % Pre-allocation
    nNeurons = length(nois);
    spks_fmk = nan(nNeurons,1001*ntrials);
    frs_fmk = nan(nNeurons,floor(1001/50)*ntrials);
    types = cell(nNeurons,1);
    regions = cell(nNeurons,1);
    channels = nan(nNeurons,1);
    clusters = nan(nNeurons,1); 
    gains = zeros(nNeurons,1);

    % Loading each neuron's data and saving some data about it
    for ii = 1:length(nois)
        load([fpath,filenames(ii+nois(1)-1).name],'type',...
            'BaselineAlignedSpikes','StopAlignedSpikes',...
            'region','channel','cluster');
        if strcmp(type,'nochange') || strcmp(type,'bad')
            gains(ii) = nan;
        else
            load([fpath,filenames(ii+nois(1)-1).name],'gain')
            gains(ii) = gain;
        end
%         smdat = smoothdata(BaselineAlignedSpikes,2,'gaussian',400);
        smdat = smoothdata(StopAlignedSpikes(:,1:1001),2,'gaussian',400);

        
        % Normalize by dividing by sqrt(max(firingrate))
        smdat = smdat / sqrt(max(smdat,[],'all'));
        
        % Normalize by taking z-score
%         smdat = zscore(smdat,[],'all');
        
%         frdat = firingRate(BaselineAlignedSpikes);
        spks_fmk(ii,:) = reshape(smdat(1:ntrials,:)',[],1);
%         frs_fmk(ii,:) = reshape(frdat(1:ntrials,:)',[],1);
        types{ii} = type;
        regions{ii} = region;
        channels(ii) = channel;
        clusters(ii) = cluster;
    end

    nCk = nchoosek(1:length(types),2);

%     spks_sm = smoothdata(spks_fmk,2,'gaussian',400);
    
    % Optional decimation so data isn't so correlated
    spks_fmk = spks_fmk(:,1:50:end);
    B = spks_fmk;

    nc_fmk = nan(length(nCk),1);
    typematch = nan(length(nCk),1);
    type1 = cell(length(nCk),1);
    type2 = cell(length(nCk),1);
    reg1 = cell(length(nCk),1);
    reg2 = cell(length(nCk),1);
    gain1 = cell(length(nCk),1);
    gain2 = cell(length(nCk),1);

    disp(length(nCk))
   
    tic
    for ii = 1:length(nCk)

        
        % First check to see if they're the same channel
        if channels(nCk(ii,1)) == channels(nCk(ii,2))
            nc_fmk(ii) = nan;
            typematch(ii) = nan;
            type1{ii} = nan;
            reg1{ii} = nan;
            reg2{ii} = nan;
            gain1{ii} = nan;
            gain2{ii} = nan;
            continue
        end
        
        % Get 2x2 matrix of correlation coefficients
        
        % Correlations of smoothed data --make sure savepath is /smoothed
        corr_temp = corr([spks_fmk(nCk(ii,1),:)' spks_fmk(nCk(ii,2),:)'],...
            'Type','Spearman');

        % Correlations of binned data --make sure savepath is /binned
%         corr_temp = corr([frs_fmk(nCk(ii,1),:)' frs_fmk(nCk(ii,2),:)'],...
%             'Type','Spearman');

        nc_fmk(ii) = corr_temp(1,2); % Get correlation of interest

        if strcmp(types{nCk(ii,1)},'inc') & strcmp(types{nCk(ii,2)},'inc')
            typematch(ii) = 1;
        elseif strcmp(types{nCk(ii,1)},'dec') & strcmp(types{nCk(ii,2)},'dec')
            typematch(ii) = 2;
        end
        if mod(ii,1000)==0,toc,end   

        type1{ii} = types{nCk(ii,1)};
        type2{ii} = types{nCk(ii,2)};

        reg1{ii} = regions{nCk(ii,1)};
        reg2{ii} = regions{nCk(ii,2)};

        gain1{ii} = gains(nCk(ii,1));
        gain2{ii} = gains(nCk(ii,2));
    end


    nc_stats = zeros(3,3);

    for ii=0:2
        if ii == 0
            x = nc_fmk(~isnan(nc_fmk));
        else
            x = nc_fmk(~isnan(nc_fmk)&typematch==ii);
        end

        SEM = std(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        CI = mean(x) + ts*SEM;

        nc_stats(ii+1,1) = mean(x);
        nc_stats(ii+1,2) = CI(1);
        nc_stats(ii+1,3) = CI(2);
    end


    X = categorical({'Inc-Dec','Inc-Inc','Dec-Dec'});
    X = reordercats(X,{'Inc-Dec','Inc-Inc','Dec-Dec'});

    Y = nc_stats(:,1);
    % Make a bar plot
    figure
    bar(X,Y)

    hold on

    errorbar([1 2 3],[nc_stats(:,1)],[nc_stats(:,3)-nc_stats(:,1)],'o')

    pwCorrs = table(nc_fmk,typematch,type1,type2,reg1,reg2,gain1,gain2);
    pwCorrs.subject = ones(size(nc_fmk))*subj;
    pwCorrs.session = ones(size(nc_fmk))*session;

    writetable(pwCorrs,['./pwCorr/smoothedCorr_chanFix/sub',num2str(subj),'ses',num2str(session),'.csv'])
end


%%

nettype = 'flux';
Ksyn2use = 0.5;

fpath = ['./data/' nettype '/Ksyn_' num2str(Ksyn2use) '/'];

files = dir([fpath,'*.mat']);

load([fpath files(1).name],'spks','time','norms')

for jj = 1:length(files)
    
    
    % Load data from first neuron in that session to pre-allocate info
    load([fpath,files(jj).name],'spks','time','type','network');
    types = type;

    % Figures path
    figures_path = 'figures/scratch';

    % Pre-allocation
    nNeurons = size(spks,2);
    ntrials = sum(~isnan(spks(:,1,1)));
    spks_fmk = nan(nNeurons,501*ntrials);
    

    % Loading each neuron's data and saving some data about it
    for ii = 1:nNeurons
                
        BaselineAlignedSpikes = squeeze(spks(:,ii,1:501));
        BaselineAlignedSpikes = BaselineAlignedSpikes(~isnan(BaselineAlignedSpikes(:,1)),:);
        
        smdat = smoothdata(BaselineAlignedSpikes,2,'gaussian',400);
        
        % Normalize by dividing by sqrt(max(firingrate))
        smdat = smdat / sqrt(max(smdat,[],'all'));
        
        % Normalize by taking z-score
%         smdat = zscore(smdat,[],'all');
        
        spks_fmk(ii,:) = reshape(smdat(1:ntrials,:)',[],1);

    end

    nCk = nchoosek(1:length(types),2);

%     spks_sm = smoothdata(spks_fmk,2,'gaussian',400);
    
    % Optional decimation so data isn't so correlated
    spks_fmk = spks_fmk(:,1:50:end);
    B = spks_fmk;

    nc_fmk = nan(length(nCk),1);
    typematch = nan(length(nCk),1);
    type1 = cell(length(nCk),1);
    type2 = cell(length(nCk),1);
    


    tic
    for ii = 1:length(nCk)

        
        % First check to see if they're the same channel
%         if channels(nCk(ii,1)) == channels(nCk(ii,2))
%             nc_fmk(ii) = nan;
%             typematch(ii) = nan;
%             type1{ii} = nan;
%             reg1{ii} = nan;
%             reg2{ii} = nan;
%             gain1{ii} = nan;
%             gain2{ii} = nan;
%             continue
%         end
        
        % Get 2x2 matrix of correlation coefficients
        
        % Correlations of smoothed data --make sure savepath is /smoothed
        corr_temp = corr([spks_fmk(nCk(ii,1),:)' spks_fmk(nCk(ii,2),:)'],...
            'Type','Spearman');

        % Correlations of binned data --make sure savepath is /binned
%         corr_temp = corr([frs_fmk(nCk(ii,1),:)' frs_fmk(nCk(ii,2),:)'],...
%             'Type','Spearman');

        nc_fmk(ii) = corr_temp(1,2); % Get correlation of interest

        if strcmp(types{nCk(ii,1)},'inc') & strcmp(types{nCk(ii,2)},'inc')
            typematch(ii) = 1;
        elseif strcmp(types{nCk(ii,1)},'dec') & strcmp(types{nCk(ii,2)},'dec')
            typematch(ii) = 2;
        end
        if mod(ii,1000)==0,toc,end   

        type1{ii} = types{nCk(ii,1)};
        type2{ii} = types{nCk(ii,2)};
        
    end


    nc_stats = zeros(3,3);

    for ii=0:2
        if ii == 0
            x = nc_fmk(~isnan(nc_fmk));
        else
            x = nc_fmk(~isnan(nc_fmk)&typematch==ii);
        end

        SEM = std(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        CI = mean(x) + ts*SEM;

        nc_stats(ii+1,1) = mean(x);
        nc_stats(ii+1,2) = CI(1);
        nc_stats(ii+1,3) = CI(2);
    end


    X = categorical({'Inc-Dec','Inc-Inc','Dec-Dec'});
    X = reordercats(X,{'Inc-Dec','Inc-Inc','Dec-Dec'});

    Y = nc_stats(:,1);
    % Make a bar plot
    figure
    bar(X,Y)

    hold on

    errorbar([1 2 3],[nc_stats(:,1)],[nc_stats(:,3)-nc_stats(:,1)],'o')

    pwCorrs = table(nc_fmk,typematch,type1,type2);
    pwCorrs.network = ones(size(nc_fmk))*network;
    

    writetable(pwCorrs,['./pwCorr/' nettype '/smoothedCorr_norm/net',num2str(network),'.csv'])
end