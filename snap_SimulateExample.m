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

% saveas(gcf,['./figures/examples/Net_',num2str(connMat_rng),'_Ksyn_',num2str(Ksyn), '_Ree_',num2str(Ree),'.svg'])
% saveas(gcf,['./figures/examples/Net_',num2str(connMat_rng),'_Ksyn_',num2str(Ksyn), '_Ree_',num2str(Ree),'.png'])



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