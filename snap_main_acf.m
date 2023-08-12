%% Main script for spontaneous fluctuations/ignitions simulations

% Main script will call the rsNet_sims function multiple times and run it
% with different parameter combinations.

% Things of interest are Ree 1:0.1:4, including sra and stp, N (400 or
% 5000), Nclusters (4 for 400, 20 and 50 for 5000), noisescale
% (0.005,0.001,0.0002)
clear all

% addpath('./functions')
% addpath('./','r')
addpath('./utils')

% savepath = 'data';
savepath = '/Volumes/GDrive/snap/data';

Ree_list = 1:0.25:7;
N_list = [400];
Nclusters_list = [4]; % Going to have to do some special work here.
noisescale_list = [0.5];
% modes = [1]; % 1=singletrial, 2=multitrial

Ksyns = [0,0.25,0.5,0.75,1];

t_iter_start = 1;
t_iter_end = 5;

% Testing
% Ree_list = 2.5;
% N_list = 400;
% Nclusters_list = 10;
% noisescale_list = [0.005,0.001];
% modes = [1,2]; % 1 for singletrial, 2 for multitrial

model = [];

tic
for i_Ksyn = 1:length(Ksyns)
    for i_noise = 1:length(noisescale_list)

        % IF Ree_vary is true, select the Ree_list specifically for
        % combinations of noise levels and synapse types

        for i_N = 1:length(N_list)
            for i_Nclust = 1:length(Nclusters_list)
                for i_iter=t_iter_start:t_iter_end % 3 simulation passes for each parameter combination.
                    for i_Ree = 1:length(Ree_list)
                            
                        %% Initialize parameters in model    

                        
                        % Get parameter values
                        Ree = Ree_list(i_Ree);
                        N = N_list(i_N);
                        Nclusters = Nclusters_list(i_Nclust);

                        noisescale = noisescale_list(i_noise);
                        Ksyn = Ksyns(i_Ksyn);
                        connMat_rng = i_iter;
                        Jee = 0.2;
                        Jclust = 1.9;


                        

                        % Create the connectivity matrix
                        W = connMatrix('N',N,...
                            'Nclusters',Nclusters,...
                            'Ree',Ree,...
                            'Jee',Jee,...
                            'Jclust',Jclust,...
                            'connMat_rng',connMat_rng);

                        % Run simulation
                        
                        model.W = W;
                        model.Ree = Ree;
                        model.N = N;
                        model.Nclusters = Nclusters;

                        model.noisescale = noisescale;
                        model.Ksyn = Ksyn;
                        model.connMat_rng = connMat_rng;
                        model.Jee = Jee;
                        model.Jclust = Jclust;


                        %% Simulate 50 trials of 1.2 seconds, calculate timescales
                        Ntrials = 20;
                        tsteps = 1200;
                        
                        spk_mats = zeros(Ntrials,N,tsteps); % Spiking behavior, Ntrials X Nneurons X time
                        taus = zeros(N,1);                    % Intrinsic timescales
                        inclusion = zeros(N,5);             % Whether the fitting resulted in problems for each neuron.
                        acf = zeros(N,12);                  % Autocorrelation
                        Rs = zeros(size(spk_mats)); 
                        parfor i=1:Ntrials
                            [spk_mats(i,:,:),Rs(i,:,:)] = snap_singletrial(W,'tsteps',1200,...
                                'noisescale',noisescale,...
                                'Ksyn',Ksyn); 

%                             spk_mats(i,:,:) = t_spk_mat; % i is in first of 3 dims because iterating over trials
%                             Rs(i,:,:) = t_rs; % the exp decayt 
                        end
                        
                        FRs = zeros(size(spk_mats,1),size(spk_mats,2),size(spk_mats,3)/50); 
                        parfor i=1:size(FRs,1)
                            FRs(i,:,:) = firingRate(squeeze(spk_mats(i,:,:)));
                        %     waitbar(i/size(FRs,1))
                        end
                        
                        % Calculate timescales for individual neurons
                        parfor i=1:N
                            [taus(i), acf(i,:), inclusion(i,:)] = calculateTimescale('spk_mat',spk_mats(:,i,201:end)); % i is in second of 3 dims because iterating over neurons
                        end
                        
                        % Calculate timescales for the cluster-averaged
                        % firing rates
                        clust_size = 0.8*N / Nclusters;
                        cluster_FRs = squeeze(mean(reshape(FRs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));
%                         PSPs = zeros(size(Rs));
%                         for i=1:Ntrials
%                             PSPs(i,:,:) = squeeze(W * squeeze(Rs(i,:,:)));
%                         end
%                         cluster_PSPs = squeeze(mean(reshape(PSPs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));
                        
                        taus_clust = zeros(Nclusters,1);                    % Intrinsic timescales
                        inclusion_clust = zeros(Nclusters,5);             % Whether the fitting resulted in problems for each neuron.
                        acf_clust = zeros(Nclusters,12);  
                        
                        parfor i=1:Nclusters
                            [taus_clust(i), acf_clust(i,:), inclusion_clust(i,:)] = calculateTimescale('FR',squeeze(cluster_FRs(:,i,5:end))); % i is in second of 3 dims because iterating over neurons
                        end
                        
                        model.Ntrials = Ntrials;
                        model.tsteps_short = tsteps;
                        model.taus_xcorr = taus;
                        model.acf_xcorr = acf;
                        model.inclusion_xcorr = inclusion;
                        
                        model.taus_clust_xcorr = taus_clust;
                        model.acf_clust_xcorr = acf_clust;
                        model.inclusion_clust_xcorr = inclusion_clust;
                        
                        
                        
                        %% Simulate one trial for 20 seconds
                        Ntrials = 1;
                        tsteps = 30000;
                        
                        spk_mats = zeros(Ntrials,N,tsteps); % Spiking behavior, Ntrials X Nneurons X time
                        taus = zeros(N,1);                    % Intrinsic timescales
                        inclusion = zeros(N,5);             % Whether the fitting resulted in problems for each neuron.
                        acf = zeros(N,12);                  % Autocorrelation
                        Rs = zeros(size(spk_mats)); 
                        parfor i=1:Ntrials
                            [spk_mats(i,:,:), Rs(i,:,:)] = snap_singletrial(model.W,'tsteps',tsteps,...
                                'noisescale',model.noisescale,...
                                'Ksyn',Ksyn); 

%                             spk_mats(i,:,:) = t_spk_mat; % i is in first of 3 dims because iterating over trials
%                             Rs(i,:,:) = t_rs; % the exp decayt 
                        end
                        
                        FRs = zeros(size(spk_mats,1),size(spk_mats,2),size(spk_mats,3)/50); 
                        for i=1:size(FRs,1)
                            FRs(i,:,:) = firingRate(squeeze(spk_mats(i,:,:)));
                        %     waitbar(i/size(FRs,1))
                        end
                        
                        % Calculate timescales for individual neurons
                        parfor i=1:N
                            [taus(i), acf(i,:), inclusion(i,:)] = calculateTimescaleACF('spk_mat',spk_mats(:,i,201:end)); % i is in second of 3 dims because iterating over neurons
                        end
                        
                        % Calculate timescales for the cluster-averaged
                        % firing rates
                        cluster_FRs = squeeze(mean(reshape(FRs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));
%                         PSPs = zeros(size(Rs));
%                         for i=1:Ntrials
%                             PSPs(i,:,:) = squeeze(W * squeeze(Rs(i,:,:)));
%                         end
%                         cluster_PSPs = squeeze(mean(reshape(PSPs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));
                        
                        taus_clust = zeros(Nclusters,1);                    % Intrinsic timescales
                        inclusion_clust = zeros(Nclusters,5);             % Whether the fitting resulted in problems for each neuron.
                        acf_clust = zeros(Nclusters,12);  
                        
                        parfor i=1:Nclusters
                            [taus_clust(i), acf_clust(i,:), inclusion_clust(i,:)] = calculateTimescaleACF('FR',cluster_FRs(i,5:end)); % i is in second of 3 dims because iterating over neurons
                        end
                        
                        model.tsteps_long = Ntrials;
                        model.taus_acf = taus;
                        model.acf_acf = acf;
                        model.inclusion_acf = inclusion;
                        
                        model.taus_clust_xcorr = taus_clust;
                        model.acf_clust_xcorr = acf_clust;
                        model.inclusion_clust_xcorr = inclusion_clust;
                        
                        
                        t_fpath = [savepath,'/noisescale_',num2str(noisescale),'/Ksyn_',num2str(Ksyn),'/'];
                        if ~exist(t_fpath,'dir')
                            mkdir(t_fpath)
                        end
                        
                        t_fname = sprintf([t_fpath,'Ree_%0.2f_rng_%d_%s.mat'],Ree,connMat_rng,datestr(now,'mm-dd-yyyy HH-MM'));
                        
                        disp(['saving to: ',t_fname])
                        
                        save(t_fname,'-struct','model')
                        
                        
                        
                        


                    end
                end
            end
        end
    end
end


toc
