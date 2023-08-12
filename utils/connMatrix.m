function [W, Wadj] = connMatrix(varargin)

%% Create a connectivity matrix for an E-I network

% By default, 80% excitatory, 20% inhibitory with balanced exciation and
% inhibition (i.e. inhibition is 4X excitation). 
% If varargin = 0, will just create a balanced asynchronous E-I matrix

% Input:

% N:            Number of neurons (default 5000)
% Pconn:        Probability of connection between any two neurons (default
% 0.2)
% cluster:      Boolean whether to cluster neurons
% Nclusters:    Number of clusters (default 500 for 80 neurons/cluster)
% cluster_inh:  Boolean whether inhibitory neurons are cluster-specific.
% Jee, Jei, Jie, Jii:   Synaptic weights for excitatory, inhibitory
% connectivity (subscripts are to, from), defaults already given
% Jclust:       Synaptic multiplier for cluster strengths
% Ree:          How heavily to bias activity towards clusters
% Spatial:      Whether to do spatially-based connectivity

% Output:

% W:            Weight matrix (0 or 1 entries)
% Wadj:         Adjacency matrix (0 or 1 entries)
% exc_mask:     Vector mask for excitatory neurons
% inh_mask:     Vector mask for inhibitory neurons
% cluster_mask: Matrix of size [N_neurons X N_clusters] where column i is a
% mask to return only the neurons of cluster i.



%% Parsing optional input

% Default param values
N = 5000;           % Number of neurons
perc_exc = 0.8;     % Percent excitatory

cluster = true;    % Whether to cluster
Nclusters = N*perc_exc/100;    % How many clusters should we use
Jee = .36;         % Default ee connectivity weights
Jei = -.675*2.5;
Jie = .21*2.5;
Jii = -0.855*2.5;
Jclust = 1.9;       % Multiplier for within-cluster synaptic weights
Pee = 0.2;          % Probability of ee connections
Pei = 0.2;
Pie = 0.2;
Pii = 0.2;

Ree = 2.5;          % ratio of in_cluster vs out_cluster connection prob, p_in/p_out.
% Note: obtain p_in and p_out using this plus the constraint that weighted average
% of p_in and p_out should be equal to p_conn.

connMat_rng = nan;


% Read in varargin
p = inputParser;
addParameter(p,'N',N,@isnumeric);
addParameter(p,'perc_exc',perc_exc,@isnumeric);
addParameter(p,'cluster',cluster,@islogical);
addParameter(p,'Nclusters',Nclusters,@isnumeric);
addParameter(p,'Jee',Jee,@isnumeric);
addParameter(p,'Jei',Jei,@isnumeric);
addParameter(p,'Jie',Jie,@isnumeric);
addParameter(p,'Jii',Jii,@isnumeric);
addParameter(p,'Jclust',Jclust,@isnumeric);
addParameter(p,'Pee',Pee,@isnumeric);
addParameter(p,'Pei',Pei,@isnumeric);
addParameter(p,'Pie',Pie,@isnumeric);
addParameter(p,'Pii',Pii,@isnumeric);
addParameter(p,'Ree',Ree,@isnumeric);
addParameter(p,'Pee_in',nan,@isnumeric);
addParameter(p,'Pee_out',nan,@isnumeric);
addParameter(p,'connMat_rng',connMat_rng);
addParameter(p,'scale',true);
addParameter(p,'syn_var',0,@isnumeric)


parse(p,varargin{:})

% Setting the variables according to default values and what is input by
% vargin

N = p.Results.N;
perc_exc = p.Results.perc_exc; 
cluster = p.Results.cluster;    % Whether to cluster
Nclusters = p.Results.Nclusters;    % How many clusters should we use
% Note: make sure Nclusters divides 0.8*N evenly!

if p.Results.scale
    Jee = p.Results.Jee * 1/sqrt(N/5000);         % Default ee connectivity weights, optimized for N=5000 and scaled accordingly
    Jei = p.Results.Jei * 1/sqrt(N/5000);
    Jie = p.Results.Jie * 1/sqrt(N/5000);
    Jii = p.Results.Jii * 1/sqrt(N/5000);
else
    Jee = p.Results.Jee;         % Default ee connectivity weights, optimized for N=5000 and scaled accordingly
    Jei = p.Results.Jei;
    Jie = p.Results.Jie;
    Jii = p.Results.Jii;
end
Jclust = p.Results.Jclust;       % Multiplier for within-cluster synaptic weights

Ree = p.Results.Ree;          % ratio of in_cluster vs out_cluster connection prob, p_in/p_out.
if cluster == false
    Ree = 1;Jclust=1;
end


% Whether we want to specify a specific RNG seed for generating matrix 
connMat_rng = p.Results.connMat_rng; 
if ~isnan(connMat_rng), rng(connMat_rng),end


%% Setting up connectivity matrix

% Major considerations:
% Are we doing clustered connectivity?
% Are we doing cluster-specific inhibition or general inhibition?




% EE connectivity:

Ne = round(perc_exc * N); % Number of excitatory neurons
Ni = round((1-perc_exc) * N); % Number of inhibitory neurons

Wee = zeros(Ne); % Pre-allocate the matrix.

cluster_size = Ne / Nclusters; % Calculate cluster size

% Calculate Pee_in and Pee_out

if isnan(p.Results.Pee_in)
    A = Nclusters * cluster_size^2;
    B = Nclusters * (Nclusters - 1) * cluster_size^2;
    Pee_in = Ne^2 * Pee * Ree / (B + A*Ree);
    Pee_out = Pee_in * (1/Ree);
else
    Pee_in = p.Results.Pee_in;
    Pee_out = p.Results.Pee_out;
    if isnan(Pee_out)
        disp(['No value set for Pee_out, using Ree=',num2str(Ree),' to set'])
        Pee_out = Pee_in * (1/Ree);
    end       
end


for i=1:Nclusters % iterate over rows
    for j=1:Nclusters % iterate over columns
        
        % Indexing: 1 + (i-1)*cluster_size iterates over first index for
        % cluster (1, 101, 201 for csize 100), i*cluster_size iterates over
        % second index for cluster (100, 200, 300 for csize 100). Same for
        % column index j.
        
        % Condition on whether or not i=j (i.e. are we in the same cluster)
        if i==j 
            
            % Condition on whether we are doing clustered connectivity
            if cluster
                Wee(1+(i-1)*cluster_size:i*cluster_size,...
                1+(j-1)*cluster_size:j*cluster_size) = ...
                    (rand(cluster_size) < Pee_in)*Jee*Jclust; 
            
            
            else 
                Wee(1+(i-1)*cluster_size:i*cluster_size,...
                1+(j-1)*cluster_size:j*cluster_size) = ...
                    (rand(cluster_size) < Pee_in)*Jee;
            end
        else
            Wee(1+(i-1)*cluster_size:i*cluster_size,1+(j-1)*cluster_size:j*cluster_size) = ...
                (rand(cluster_size) < Pee_out)*Jee;
        end      
    end
end

% IE connectivity

Wie = (rand(Ni,Ne) < Pie) * Jie;

% EI connectivity
Wei = (rand(Ne,Ni) < Pei) * Jei;

% II connectivity

Wii = (rand(Ni,Ni) < Pii) * Jii;


W = [Wee,Wei;Wie,Wii];

W = W - diag(diag(W));

Wadj = W~=0;

% if we specified a random seed when generating the matrix, shuffle the rng
if ~isnan(connMat_rng),rng shuffle,end
    
end