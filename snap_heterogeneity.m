%% SNAP Project Heterogeneity Analysis

% Here we test two models of heterogeneity in how neurons 'ramp up' their
% firing rates. First is a linear increase / ramp from ~5 Hz to ~20 Hz and
% then next is a stepping increase/ramp from the same but with a
% normally-distributed ramping 'onset' time. 



%% Linear Ramping Model -- Increasing

nSims = 100; % Number of total simulations
nTrials = 50; % number of trials per simulated neuron
time = -2.95:0.05:0.;


norms_inc_linear = zeros(nSims,1);

cv2_inc_linear = zeros(nSims,1);

ff_inc_linear = zeros(nSims,1);


for jj = 1:nSims
    this_fr_start = abs(5 + randn);

    this_fr_end = 20 + randn;

    frs = linspace(this_fr_start,this_fr_end,length(time));

    spks = zeros(nTrials,length(time)*50);



    for kk =1:length(time)
        spks(:,1 + (kk-1)*50:kk*50) = rand(size(spks(:,1+(kk-1)*50:kk*50))) < frs(kk)/1000;
    end


    spks_sm = smoothdata(spks,2,'gaussian',400)*1000;
    
    spks_sm = spks_sm / max(spks_sm,[],'all');

    spks_sm_mean = nanmean(spks_sm,1);
    
    nCk = nchoosek(1:nTrials,2);

    norms_this = 0;
    for kk=1:length(nCk)
%         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
        norms_this = norms_this + norm(spks_sm(nCk(kk,1))-spks_sm(nCk(kk,2)));

    end
    
    norms_inc_linear(jj) = norms_this / length(nCk);
    
    A=[];
    for ii = 1:size(spks,1)
        A = [A calculateCV2(spks(ii,:))];
    end
    cv2_inc_linear(jj) = mean(A(A ~= 0));
    
    ff_inc_linear(jj) = calculateFanoFactor(spks(:,time < -1));
end



% Plot the last one
plot(linspace(-3,0,size(spks_sm,2)),spks_sm)
hold on
plot(linspace(-3,0,size(spks_sm,2)),spks_sm_mean,'k','LineWidth',2)

%% Ramping Decreasing Model

nSims = 100; % Number of total simulations
nTrials = 50; % number of trials per simulated neuron
time = -2.95:0.05:0.;


norms_dec_linear = zeros(nSims,1);
cv2_dec_linear = zeros(nSims,1);

ff_dec_linear = zeros(nSims,1);

for jj = 1:nSims
    this_fr_start = abs(20 + randn);

    this_fr_end = abs(5 + randn);

    frs = linspace(this_fr_start,this_fr_end,length(time));

    spks = zeros(nTrials,length(time)*50);



    for kk =1:length(time)
        spks(:,1 + (kk-1)*50:kk*50) = rand(size(spks(:,1+(kk-1)*50:kk*50))) < frs(kk)/1000;
    end


    spks_sm = smoothdata(spks,2,'gaussian',400)*1000;
    
    spks_sm = spks_sm / max(spks_sm,[],'all');


    spks_sm_mean = nanmean(spks_sm,1);

    nCk = nchoosek(1:nTrials,2);

    norms_this = 0;
    for kk=1:length(nCk)
%         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
        norms_this = norms_this + norm(spks_sm(nCk(kk,1))-spks_sm(nCk(kk,2)));

    end
    
    norms_dec_linear(jj) = norms_this / length(nCk);
    
    A=[];
    for ii = 1:size(spks,1)
        A = [A calculateCV2(spks(ii,:))];
    end
    cv2_dec_linear(jj) = mean(A(A ~= 0));
    
    ff_dec_linear(jj) = calculateFanoFactor(spks(:,time < -1));

end

% Plot the last one
figure
plot(linspace(-3,0,size(spks_sm,2)),spks_sm)
hold on
plot(linspace(-3,0,size(spks_sm,2)),spks_sm_mean,'k','LineWidth',2)
%% Stepping Model Increasing
nSims = 100; % Number of total simulations
nTrials = 50; % number of trials per simulated neuron
time = -2.95:0.05:0.;


norms_inc_step = zeros(nSims,1);

cv2_inc_step = zeros(nSims,1);

ff_inc_step = zeros(nSims,1);

for jj = 1:nSims
    this_fr_start = abs(5 + randn);

    this_fr_end = 20 + randn;
    
    this_ramp_mean = -1.5;
    this_ramp_var = 0.5;

    frs = ones(length(time),1);
    

    spks = zeros(nTrials,length(time)*50);



    for ll = 1:nTrials
        
        this_ramp_onset = this_ramp_mean + randn * this_ramp_var;
        for kk =1:length(time)
            if time(kk) < this_ramp_onset
                spks(ll,1 + (kk-1)*50:kk*50) = rand(size(spks(ll,1+(kk-1)*50:kk*50))) < this_fr_start/1000;
            else
                spks(ll,1 + (kk-1)*50:kk*50) = rand(size(spks(ll,1+(kk-1)*50:kk*50))) < this_fr_end/1000;
            end
        end
    end


    spks_sm = smoothdata(spks,2,'gaussian',400)*1000;
    
    spks_sm = spks_sm / max(spks_sm,[],'all');

    spks_sm_mean = nanmean(spks_sm,1);

    nCk = nchoosek(1:nTrials,2);

    norms_this = 0;
    for kk=1:length(nCk)
%         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
        norms_this = norms_this + norm(spks_sm(nCk(kk,1))-spks_sm(nCk(kk,2)));

    end
    
    norms_inc_step(jj) = norms_this / length(nCk);
    
    A=[];
    for ii = 1:size(spks,1)
        A = [A calculateCV2(spks(ii,:))];
    end
    cv2_inc_step(jj) = mean(A(A ~= 0));
    
    ff_inc_step(jj) = calculateFanoFactor(spks(:,time < -1));

end
% Plot the last one
figure
plot(linspace(-3,0,size(spks_sm,2)),spks_sm)
hold on
plot(linspace(-3,0,size(spks_sm,2)),spks_sm_mean,'k','LineWidth',2)

%% Stepping Model Decreasing
nSims = 100; % Number of total simulations
nTrials = 50; % number of trials per simulated neuron
time = -2.95:0.05:0.;


norms_dec_step = zeros(nSims,1);
cv2_dec_step = zeros(nSims,1);

ff_dec_step = zeros(nSims,1);

for jj = 1:nSims
    this_fr_start = 20 + randn;

    this_fr_end = abs(5 + randn);
    
    this_ramp_mean = -1.5;
    this_ramp_var = 0.5;

    frs = ones(length(time),1);
    

    spks = zeros(nTrials,length(time)*50);



    for ll = 1:nTrials
        
        this_ramp_onset = this_ramp_mean + randn * this_ramp_var;
        for kk =1:length(time)
            if time(kk) < this_ramp_onset
                spks(ll,1 + (kk-1)*50:kk*50) = rand(size(spks(ll,1+(kk-1)*50:kk*50))) < this_fr_start/1000;
            else
                spks(ll,1 + (kk-1)*50:kk*50) = rand(size(spks(ll,1+(kk-1)*50:kk*50))) < this_fr_end/1000;
            end
        end
    end


    spks_sm = smoothdata(spks,2,'gaussian',400)*1000;
    spks_sm = spks_sm / max(spks_sm,[],'all');


    spks_sm_mean = nanmean(spks_sm,1);

    nCk = nchoosek(1:nTrials,2);

    norms_this = 0;
    for kk=1:length(nCk)
%         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
        norms_this = norms_this + norm(spks_sm(nCk(kk,1))-spks_sm(nCk(kk,2)));

    end
    
    norms_dec_step(jj) = norms_this / length(nCk);
    
    A=[];
    for ii = 1:size(spks,1)
        A = [A calculateCV2(spks(ii,:))];
    end
    cv2_dec_step(jj) = mean(A(A ~= 0));
    
    ff_dec_step(jj) = calculateFanoFactor(spks(:,time < -1));

end
% Plot the last one
figure
plot(linspace(-3,0,size(spks_sm,2)),spks_sm)
hold on
plot(linspace(-3,0,size(spks_sm,2)),spks_sm_mean,'k','LineWidth',2)

%% Load in the KSYN values


nettype = 'flux'; % flux or ignition

fpath = ['./data/' nettype '/Ksyn_' num2str(0.5) '/'];
files = dir([fpath,'*.mat']);


norms_inc_ksyn5 = zeros(length(files),1);
norms_dec_ksyn5 = zeros(length(files),1);

% cv2_inc_ksyn5 = cell(length(files),1);
% cv2_dec_ksyn5 = cell(length(files),1);

ff_inc_ksyn5 = cell(length(files),1);
ff_dec_ksyn5 = cell(length(files),1);

tic
for ii = 1:length(files)
    
    % Loading for pre-allocation
    load([fpath files(ii).name],'norms','type','spks')
%     load([fpath files(ii).name],'type','spks')

    norms_inc_ksyn5(ii) = nanmedian(norms(strcmp(type,'inc')));
    norms_dec_ksyn5(ii) = nanmedian(norms(strcmp(type,'dec')));
    
    trials2use = ~isnan(spks(:,1,1));

    
    cv2_inc_this = nan(320,1);
    cv2_dec_this = nan(320,1);
    
    ff_inc_this = nan(320,1);
    ff_dec_this = nan(320,1);
    parfor jj = 1:size(spks,2)
        
        
        
        if strcmp(type{jj},'inc')
%             A = [];
%             for kk = 1:size(spks,1)
%                 A = [A calculateCV2(squeeze(spks(kk,jj,:)))];
%             end
%             cv2_inc_this(jj) = nanmean(A(A ~= 0));
            
            ff_inc_this(jj) = calculateFanoFactor(squeeze(spks(trials2use,jj,time < -1)));
        elseif strcmp(type{jj},'dec')
%             A=[];
%             for kk = 1:size(spks,1)
%                 A = [A calculateCV2(squeeze(spks(kk,jj,:)))];
%             end
%             cv2_dec_this(jj) = nanmean(A(A ~= 0));
            
            ff_dec_this(jj) = calculateFanoFactor(squeeze(spks(trials2use,jj,time < -1)));
        end
    
    end
    
%     cv2_inc_ksyn5{ii} = cv2_inc_this(~isnan(cv2_inc_this));
%     cv2_dec_ksyn5{ii} = cv2_dec_this(~isnan(cv2_dec_this));
    
    ff_inc_ksyn5{ii} = ff_inc_this(~isnan(ff_inc_this));
    ff_dec_ksyn5{ii} = ff_dec_this(~isnan(ff_dec_this));
    toc
end


fpath = ['./data/' nettype '/Ksyn_' num2str(0) '/'];
files = dir([fpath,'*.mat']);


norms_inc_ksyn0 = nan(length(files),1);
norms_dec_ksyn0 = nan(length(files),1);

% cv2_inc_ksyn0 = cell(length(files),1);
% cv2_dec_ksyn0 = cell(length(files),1);

ff_inc_ksyn0 = cell(length(files),1);
ff_dec_ksyn0 = cell(length(files),1);

for ii = 1:length(files)
    % Loading for pre-allocation
    
    try
        load([fpath files(ii).name],'norms','type','spks')
%         load([fpath files(ii).name],'type','spks')
    catch
        continue
    end
    
    trials2use = ~isnan(spks(:,1,1));
    
    norms_inc_ksyn0(ii) = nanmedian(norms(strcmp(type,'inc')));
    norms_dec_ksyn0(ii) = nanmedian(norms(strcmp(type,'dec')));
    
    cv2_inc_this = nan(320,1);
    cv2_dec_this = nan(320,1);
    
    ff_inc_this = nan(320,1);
    ff_dec_this = nan(320,1);
    parfor jj = 1:size(spks,2)
        
        if strcmp(type{jj},'inc')
%             A = [];
%             for kk = 1:size(spks,1)
%                 A = [A calculateCV2(squeeze(spks(kk,jj,:)))];
%             end
%             cv2_inc_this(jj) = nanmean(A(A ~= 0));
            ff_inc_this(jj) = calculateFanoFactor(squeeze(spks(trials2use,jj,time < -1)));

        elseif strcmp(type{jj},'dec')
%             A=[];
%             for kk = 1:size(spks,1)
%                 A = [A calculateCV2(squeeze(spks(kk,jj,:)))];
%             end
%             cv2_dec_this(jj) = nanmean(A(A ~= 0));
            ff_dec_this(jj) = calculateFanoFactor(squeeze(spks(trials2use,jj,time < -1)));

        end
    
    end
    
%     cv2_inc_ksyn0{ii} = cv2_inc_this(~isnan(cv2_inc_this));
%     cv2_dec_ksyn0{ii} = cv2_dec_this(~isnan(cv2_dec_this));
    ff_inc_ksyn5{ii} = ff_inc_this(~isnan(ff_inc_this));
    ff_dec_ksyn5{ii} = ff_dec_this(~isnan(ff_dec_this));
    toc
end


%% Real CV2 from FMK data


fpath = ('data_real/');

filenames = dir([fpath,'*.mat']);

load([fpath,filenames(1).name]);

cv2_inc_real = nan(length(filenames),1);
cv2_dec_real = nan(length(filenames),1);

ff_inc_real = nan(length(filenames),1);
inc_gain = nan(length(filenames),1);
ff_dec_real = nan(length(filenames),1);
dec_gain = nan(length(filenames),1);

for ii = 1:length(filenames)
    load([fpath,filenames(ii).name],'type','StopAlignedSpikes','gain')
    
    if ~(strcmp(type,'inc')||strcmp(type,'dec')),continue,end
    
%     A = [];
%     for kk = 1:size(StopAlignedSpikes,1)
%         A = [A calculateCV2(StopAlignedSpikes(kk,:))];
%     end
    
    if strcmp(type,'inc')
%         cv2_inc_real(ii) = mean(A(A ~= 0));
        ff_inc_real(ii) = calculateFanoFactor(StopAlignedSpikes(:,501:2500));
        inc_gain(ii) = gain;
    elseif strcmp(type,'dec')
%         cv2_dec_real(ii) = mean(A(A ~= 0));
        ff_dec_real(ii) = calculateFanoFactor(StopAlignedSpikes(:,501:2500));
        dec_gain(ii) = gain;
    end
end

cv2_inc_real = cv2_inc_real(~isnan(cv2_inc_real));
cv2_dec_real = cv2_dec_real(~isnan(cv2_dec_real));

ff_inc_real = ff_inc_real(~isnan(ff_inc_real));
ff_dec_real = ff_dec_real(~isnan(ff_dec_real));

inc_gain = inc_gain(~isnan(inc_gain));
dec_gain = dec_gain(~isnan(dec_gain));

    

%% Saving Everything


% save('./heterogeneity/l2norms.mat',...
%     'norms_inc_linear','norms_dec_linear',...
%     'norms_inc_step','norms_dec_step',...
%     'norms_inc_ksyn5','norms_dec_ksyn5',...
%     'norms_inc_ksyn0','norms_dec_ksyn0','-append')

% save('./heterogeneity/cv2.mat',...
%     'cv2_inc_linear','cv2_dec_linear',...
%     'cv2_inc_step','cv2_dec_step',...
%     'cv2_inc_ksyn5','cv2_dec_ksyn5',...
%     'cv2_inc_ksyn0','cv2_dec_ksyn0',...
%     'cv2_inc_real','cv2_dec_real')

save('./heterogeneity/ff.mat',...
    'ff_inc_linear','ff_dec_linear',...
    'ff_inc_step','ff_dec_step',...
    'ff_inc_ksyn5','ff_dec_ksyn5',...
    'ff_inc_ksyn0','ff_dec_ksyn0',...
    'ff_inc_real','ff_dec_real')
    


%% Plotting

load('./heterogeneity/cv2.mat')


%% Real vs our model


figure('position',[0 0 500 250])

tiledlayout(1,2,'TileSpacing','compact')
nexttile

% histogram(cv2_inc_real(inc_gain < nanmedian(inc_gain)),...
histogram(cv2_inc_real,...
    'FaceColor',[0.5 0.15 0.1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')
hold on
histogram(cat(1,cv2_inc_ksyn5{:}),...
    'FaceColor',[1 0.3 0.2],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')

xline(mean(cv2_inc_real),'Color',[0.5 0.15 0.1],'LineWidth',3)
xline(mean(cat(1,cv2_inc_ksyn5{:})),'Color',[1 0.3 0.2],'LineWidth',3)


yl = ylim;

text(mean(cv2_inc_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0.5 0.15 0.1],'FontSize',14)
% text(mean(cat(1,cv2_inc_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)
text(-0.08,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)

box off

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)

nexttile
% histogram(cv2_dec_real(dec_gain < nanmedian(dec_gain)),...
histogram(cv2_dec_real,...
    'FaceColor',[0 0.25 0.5],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')
hold on
histogram(cat(1,cv2_dec_ksyn5{:}),...
    'FaceColor',[0 0.5 1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')

xline(mean(cv2_dec_real),'Color',[0 0.25 0.5],'LineWidth',3)
xline(mean(cat(1,cv2_dec_ksyn5{:})),'Color',[0 0.5 1],'LineWidth',3)

yl = ylim;
text(mean(cv2_dec_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0 0.25 0.5],'FontSize',14)
% text(mean(cat(1,cv2_dec_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)
text(0,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)

ax.Background.Visible = 'off';

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)

box off

% saveas(gcf,['figures/model_preds/' nettype '_CV2_narrow.svg'])
% saveas(gcf,['figures/model_preds/' nettype '_CV2_narrow.png'])


%% Linear vs step vs real

figure('position',[0 0 800 400])

tiledlayout(1,2,'TileSpacing','compact')
nexttile

histogram(cv2_inc_real,...
    'FaceColor','black','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')
hold on
histogram(cv2_inc_linear,...
    'FaceColor','red','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')
histogram(cv2_inc_step,...
    'FaceColor','magenta','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)

nexttile

histogram(cv2_dec_real,...
    'FaceColor','black','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')
hold on
histogram(cv2_dec_linear,...
    'FaceColor','blue','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')

histogram(cv2_dec_step,...
    'FaceColor','magenta','EdgeAlpha',0,'FaceAlpha',0.5,...
    'BinEdges',0:0.05:2,'Normalization','probability')

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)


%% CV2 Changes 1 second before movement? Threshold-crossing analysis

breakpoint = -1000;

fpath = ('data_real/');

filenames = dir([fpath,'*.mat']);

load([fpath,filenames(1).name]);

cv2_inc_real_pre = nan(length(filenames),1);
cv2_dec_real_pre = nan(length(filenames),1);

cv2_inc_real_post = nan(length(filenames),1);
cv2_dec_real_post = nan(length(filenames),1);

for ii = 1:length(filenames)
    load([fpath,filenames(ii).name],'type','StopAlignedSpikes')
    
    if ~(strcmp(type,'inc')||strcmp(type,'dec')),continue,end
    
    A = [];
    B = [];
    for kk = 1:size(StopAlignedSpikes,1)
        A = [A calculateCV2(StopAlignedSpikes(kk,1:2501))];
        B = [B calculateCV2(StopAlignedSpikes(kk,2501:end))];
    end
    
    if strcmp(type,'inc')
        cv2_inc_real_pre(ii) = nanmean(A(A ~= 0));
        cv2_inc_real_post(ii) = nanmean(B(B ~= 0));
    elseif strcmp(type,'dec')
        cv2_dec_real_pre(ii) = nanmean(A(A ~= 0));
        cv2_dec_real_post(ii) = nanmean(B(B ~= 0));
    end
end

cv2_inc_real_pre = cv2_inc_real_pre(~isnan(cv2_inc_real_post));
cv2_dec_real_pre = cv2_dec_real_pre(~isnan(cv2_dec_real_post));

cv2_inc_real_post = cv2_inc_real_post(~isnan(cv2_inc_real_post));
cv2_dec_real_post = cv2_dec_real_post(~isnan(cv2_dec_real_post));

%% Pre vs post histogram
figure('position',[0 0 500 250])

tiledlayout(1,2,'TileSpacing','compact')
nexttile

histogram(cv2_inc_real_pre,...
    'FaceColor',[0.5 0.15 0.1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')
hold on
histogram(cv2_inc_real_post,...
    'FaceColor',[1 0.3 0.2],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')

xline(mean(cv2_inc_real_pre),'Color',[0.5 0.15 0.1],'LineWidth',3)
xline(mean(cv2_inc_real_post),'Color',[1 0.3 0.2],'LineWidth',3)


yl = ylim;

text(mean(cv2_inc_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0.5 0.15 0.1],'FontSize',14)
% text(mean(cat(1,cv2_inc_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)
text(-0.08,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)

box off

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)

nexttile

histogram(cv2_dec_real_pre,...
    'FaceColor',[0 0.25 0.5],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')
hold on
histogram(cat(1,cv2_dec_real_post),...
    'FaceColor',[0 0.5 1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.1:2,'Normalization','probability')

xline(mean(cv2_dec_real_pre),'Color',[0 0.25 0.5],'LineWidth',3)
xline(mean(cv2_dec_real_post),'Color',[0 0.5 1],'LineWidth',3)

yl = ylim;
% text(mean(cv2_dec_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0 0.25 0.5],'FontSize',14)
% % text(mean(cat(1,cv2_dec_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)
% text(0,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)

ax.Background.Visible = 'off';

% xticks([0 50 100])
xlabel('CV2','FontSize',16)
ylabel('Proportion','FontSize',16)

box off





%% Difference plot
figure('position',[0 0 500 250])

% tiledlayout(1,2,'TileSpacing','compact')
% nexttile
subplot(1,2,1)
for ii = 1:length(cv2_inc_real_post)
    
    plot([0 1],[cv2_inc_real_pre,cv2_inc_real_post],'Color',[1 0.3 0.2])
    hold on

end

plot([0 1],[nanmean(cv2_inc_real_pre) nanmean(cv2_inc_real_post)],'k','LineWidth',3)
[H,P,CI] = ttest(cv2_inc_real_post,cv2_inc_real_pre);
text(0.2,0.2,['p < 0.001'],'Color',[1 0.3 0.2],'FontSize',14)

title('Increasing Neurons')
xlim([-0.2 1.2]),ylim([0 1.6])
ax = gca;
axis off
ax.YAxis.Visible = 'on';
ax.XAxis.Visible = 'on';
xticks([0 1])
xticklabels({'Before -1s','After -1s'})

ylabel('CV2')

subplot(1,2,2)
for ii = 1:length(cv2_dec_real_post)
    
    plot([0 1],[cv2_dec_real_pre,cv2_dec_real_post],'Color',[0 0.5 1])
    hold on

end

plot([0 1],[nanmean(cv2_dec_real_pre) nanmean(cv2_dec_real_post)],'k','LineWidth',3)
[H,P,CI] = ttest(cv2_dec_real_post,cv2_dec_real_pre);
text(0.2,0.2,['p < 0.001'],'Color',[0 0.5 1],'FontSize',14)
title('Decreasing Neurons')
xlim([-0.2 1.2]),ylim([0 1.6])
ax = gca;
axis off
ax.YAxis.Visible = 'on';
ax.XAxis.Visible = 'on';
xticks([0 1])
xticklabels({'Before -1s','After -1s'})
% ylabel('CV2')


% saveas(gcf,['figures/model_preds/CV2_ignition.svg'])
% saveas(gcf,['figures/model_preds/CV2_ignition.png'])


%% Plotting for Fano Factors

load('./heterogeneity/ff.mat')

figure('position',[0 0 500 250])

tiledlayout(1,2,'TileSpacing','compact')
nexttile

histogram(ff_inc_real,...
    'FaceColor',[0.5 0.15 0.1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.25:10,'Normalization','probability')
hold on
histogram(cat(1,ff_inc_ksyn5{:}),...
    'FaceColor',[1 0.3 0.2],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.25:10,'Normalization','probability')

xline(mean(ff_inc_real),'Color',[0.5 0.15 0.1],'LineWidth',3)
xline(mean(cat(1,ff_inc_ksyn5{:})),'Color',[1 0.3 0.2],'LineWidth',3)


yl = ylim;

text(mean(ff_inc_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0.5 0.15 0.1],'FontSize',14)
% text(mean(cat(1,ff_inc_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)
text(-0.08,yl(2)-0.05,{'Simulated','Neurons'},'Color',[1 0.3 0.2],'FontSize',14)

box off

% xticks([0 50 100])
xlabel('Fano Factor','FontSize',16)
ylabel('Proportion','FontSize',16)

nexttile

histogram(ff_dec_real,...
    'FaceColor',[0 0.25 0.5],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.25:10,'Normalization','probability')
hold on
histogram(cat(1,ff_dec_ksyn5{:}),...
    'FaceColor',[0 0.5 1],'EdgeAlpha',0,'FaceAlpha',0.3,...
    'BinEdges',0:0.25:10,'Normalization','probability')

xline(mean(ff_dec_real),'Color',[0 0.25 0.5],'LineWidth',3)
xline(mean(cat(1,ff_dec_ksyn5{:})),'Color',[0 0.5 1],'LineWidth',3)

yl = ylim;
text(mean(ff_dec_real)+0.15,yl(2)-0.05,{'Real','Neurons'},'Color',[0 0.25 0.5],'FontSize',14)
% text(mean(cat(1,ff_dec_ksyn5{:}))-0.55,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)
text(0,yl(2)-0.05,{'Simulated','Neurons'},'Color',[0 0.5 1],'FontSize',14)

ax.Background.Visible = 'off';

% xticks([0 50 100])
xlabel('Fano Factor','FontSize',16)
ylabel('Proportion','FontSize',16)

box off

saveas(gcf,['figures/model_preds/' nettype '_FF_narrow.svg'])
saveas(gcf,['figures/model_preds/' nettype '_FF_narrow.png'])