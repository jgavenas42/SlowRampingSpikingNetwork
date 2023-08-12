%% Outline

% 1) load in the data from the 500-second simulations
%       a) rule out clusters that don't have spontaneous ignitions
%       b) new idea: decide on Ws (rng, Ree, etc.) that lead to 1/f
%       behavior before simulating the 500 second stuff.
% 2) iterate over cluster-averaged firing rates, looking for
% threshold-crossings in each cluster according to
%       a) variable threshold (30 Hz? 40 Hz? 20 Hz?)
%       b) variable re-set distance from average (w/in 3 Hz?)
% 3) find the "best" combination of threshold, re-set criteria, come up
% with a list of "threshold-crossing times"
% 4) create new data structures (ignition-aligned), and load data into
% those ignition structures
% 5) backwards-average individual neurons, clusters' FRs, and average PSPs
% according to those ignition-aligned structures.

clear all
addpath(genpath('./utils'))

ksyn_str = 'Ksyn_0.75/';
% Select our filepath
t_fpath= ['/Volumes/GDrive/snap/data_th_ex/noisescale_0.5/',ksyn_str]; % Select the filepath
% t_savepath = ['/Volumes/GDrive/snap/data_th/',ksyn_str];
t_savepath = ['data_th/',ksyn_str];

t_fid = fopen([t_fpath,'filesselected_auto.txt'],'a+'); 
t_A = fileread([t_fpath,'filesselected_auto.txt']); 
% Open and read the file where we've selected what matrices to run.
%

t_files = regexp(t_A, '\r\n|\r|\n', 'split'); % Splits

%% Run Simulations
for i=1:length(t_files)-1
% for i=1:1
    disp(t_files{i})
    
%     t_idx = strfind(t_A,t_fname);

%     if isnan('thresh_ba') %|| ~isempty(t_idx),
%         continue,
%     end
    tic
    load([t_fpath,t_files{i}])
    
    
    
                        
    model.W = W;
    model.Ree = Ree;
    model.N = N;
    model.Nclusters = Nclusters;

    model.noisescale = noisescale;
    model.Ksyn = Ksyn;
    model.connMat_rng = connMat_rng;
    model.Jee = Jee;
    model.Jclust = Jclust;


    %% Simulate 100 trials of 10 seconds
    Ntrials = 100;
    tsteps = 10000;
    
    % Pre-allocate the spike and synaptic potential variables
    spk_mats = zeros(Ntrials,N,tsteps); % Spiking behavior, Ntrials X Nneurons X time
    Rs = zeros(size(spk_mats)); 
    
    % Set the first however many trials to be from the example simulations
    spk_mats(1:ex_Ntrials,:,:) = ex_spks;
    Rs(1:ex_Ntrials,:,:) = ex_Rs;
    
    
    parfor j=(ex_Ntrials+1):Ntrials
        [spk_mats(j,:,:), Rs(j,:,:)] = snap_singletrial(W,'tsteps',tsteps,...
            'noisescale',noisescale,...
            'Ksyn',Ksyn); 
    end
    
    % Pre-allocate firing rate variable
    FRs = zeros(size(spk_mats,1),size(spk_mats,2),size(spk_mats,3)/50);
    FRs(1:ex_Ntrials,:,:) = ex_FRs; %similarly, make it from the example one
    for j=(ex_Ntrials+1):size(FRs,1)
        FRs(j,:,:) = firingRate(squeeze(spk_mats(j,:,:)));
    %     waitbar(i/size(FRs,1))
    end



%     % Calculate timescales for the cluster-averaged
%     % firing rates
    clust_size = 0.8*N / Nclusters;
    cluster_FRs = squeeze(mean(reshape(FRs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));
    PSPs = zeros(size(Rs));
    for j=1:Ntrials
        PSPs(j,:,:) = squeeze(W * squeeze(Rs(j,:,:)));
    end
    cluster_PSPs = squeeze(mean(reshape(PSPs(:,1:0.8*N,:),Ntrials,clust_size,Nclusters,[]),2));


    model.Ntrials = Ntrials;
    model.tsteps = tsteps;


    model.spks = cell(N,1);
    for j=1:N
        model.spks{j} = sparse(squeeze(spk_mats(:,j,:)));
    end
    
    model.Rs = Rs;
    model.FRs = FRs;
    model.cluster_FRs = cluster_FRs;
    model.clust_size = clust_size;
%     model.PSPs = PSPs;
%     model.cluster_PSPs = cluster_PSPs;
    
    
    if ~exist(t_savepath,'dir')
        mkdir(t_savepath)
    end

    t_fname = sprintf([t_savepath,'Ree_%0.2f_rng_%d_%s.mat'],Ree,connMat_rng,datestr(now,'mm-dd-yyyy HH-MM'));

    disp(['saving to: ',t_fname])

    save(t_fname,'-struct','model','-v7.3')
    
    
    
    toc
    
end
fclose(t_fid);

%% Plotting/Inspecting Simulated Networks

