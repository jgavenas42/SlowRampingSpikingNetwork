function [spk_mat, rs, REC] = snap_singletrial(W,varargin)

% Simulate a single continuous trial of a LIF network

%% Parsing input

% Default input values
Iapp = 300; % Applied current in pA. Should be scalar, vector size Nx1, vector size (1xtsteps) or matrix size (Nxtsteps)
tsteps = 1200;
noisescale = 1; % standard deviation of membrane noise (resting -0.07)
noiseduration = 10; % how many timesteps during which to update the noise

curr_type = 'same'; % Type of current to administer. Can be 'same','supra', or 'sub'

% Spike rate adaptation
% Constants taken from Treves (1993) Mean Field Analysis of Neuronal Spike
% Dynamics. General idea is to include a potassium-driven current that
% limits spiking.
sra = false; % Spike rate adaptation
V_k = -0.085; % Reversal potential at -85 mV
tau_k = 0.080; % 80 ms decay time-constant (changed)
delta_g_k = 7e-10; % Increment for potassium conductance upon a spike originally 9e-9

% Short term plasticity values
stp = false; % short-term plasticity in the connectivities
stp_oldmode = false;
stp_tau = 0.02; % Decay time for STP
stp_frac = 0.1; % percentage each spike decreases synaptic efficiency

tau_rise = .002;   % Synaptic rise time
tau_fall = .02;  % Synaptic fall time

% Alternative to specifying the tau_rise and tau_fall, specifying Ksyn will
% make it so that the function maintains activity of an average of fast and
% slow synapses, and Ksyn determines how they are averaged.
Ksyn = 0.0; 

p = inputParser;
addParameter(p,'Iapp',Iapp,@isnumeric);
addParameter(p,'tsteps',tsteps,@isnumeric);
addParameter(p,'stp',stp,@islogical);
addParameter(p,'stp_frac',stp_frac,@isnumeric);
addParameter(p,'stp_tau',stp_tau,@isnumeric);
addParameter(p,'sra',sra,@islogical);
addParameter(p,'tau_k',tau_k,@isnumeric);
addParameter(p,'delta_g_k',delta_g_k,@isnumeric);

addParameter(p,'noisescale',noisescale,@isnumeric);
addParameter(p,'noiseduration',noiseduration,@isnumeric);

addParameter(p,'stp_oldmode',stp_oldmode,@islogical);
addParameter(p,'curr_type',curr_type)
addParameter(p,'init_state',nan,@isnumeric)
addParameter(p,'tau_rise',tau_rise,@isnumeric)
addParameter(p,'tau_fall',tau_fall,@isnumeric)
addParameter(p,'Ksyn',Ksyn,@isnumeric)
% addParameter(p,'stim_th',nan,@isnumeric)

parse(p,varargin{:})

Iapp = p.Results.Iapp;
tsteps = p.Results.tsteps;
noisescale = p.Results.noisescale; % standard deviation of membrane noise (resting -0.07)
noiseduration = p.Results.noiseduration;

sra = p.Results.sra; % Spike rate adaptation
V_k = -0.085; % Reversal potential at -85 mV
tau_k = 0.080; % 80 ms decay time-constant (changed)
delta_g_k = 2e-10; % Increment for potassium conductance upon a spike.
stp = p.Results.stp; % short-term plasticity in the connectivities

stp_oldmode = p.Results.stp_oldmode;
stp_tau = p.Results.stp_tau; % Decay time for STP
stp_frac = p.Results.stp_frac; % percentage each spike decreases synaptic efficiency

if ~isnan(p.Results.init_state),rng(p.Results.init_state),end

Ksyn = p.Results.Ksyn;

N = size(W,1);
dt = 0.001;


% Preallocation of memory
spk_mat = zeros(size(W,1),tsteps);


Vm = -0.070 * ones(size(W,1),1); % For keeping track of membrane potential
Vm = Vm + 0.02 * rand(size(Vm)); % randomly initialize membrane potential
spks = zeros(size(W,1),1); % For keeping track of spiking.
ref = zeros(size(Vm,1),1); % For keeping track of refractory periods
ref_period = 2;

% Set up the synaptic stuff: for now we will use a double exponential
% filter. 


if ~isnan(Ksyn) % I.e. if we specify Ksyn, then we hardcode the synapses and keep track of both fast and slow synapses.
    tau_rise_fast = 0.001;
    tau_fall_fast = 0.005;
    
    tau_rise_slow = 0.020;
    tau_fall_slow = 0.100;
    
    si_fast = zeros(size(W,1),1);
    ri_fast = zeros(size(W,1),1);
    si_slow = zeros(size(W,1),1);
    ri_slow = zeros(size(W,1),1);  
else
    si = zeros(size(W,1),1);
    ri = zeros(size(W,1),1);
    tau_rise = p.Results.tau_rise;   % Synaptic rise time
    tau_fall = p.Results.tau_fall;  % Synaptic fall time
end
%% Short-term plasticity

if stp
    W_stp = zeros(size(W,1),1);
    W_stp_idx = zeros(tsteps,1);
end

if sra
    g_k = zeros(size(W,1),1);
end
    

% example_neuron = zeros(4,tsteps); % Plotting membrane voltage of a single neuron.
% example_neuron_idx = 1;

syn_curr_multiplier = 1; % Need to multiply the synaptic currents by some capacitance-like value in order to keep them on the same level.

% bias_current = Iapp/10 * [ones(cluster_size*3,1);zeros(5000-cluster_size*3,1)];
% Iapp = Iapp + Iapp * [rand([4000,1])*0.2;rand([1000,1])*0.05];
if strcmp(curr_type,'supra')
    Iapp = Iapp + Iapp * [rand([5000,1])*0.1];
elseif strcmp(curr_type,'sub')
    Iapp = Iapp - Iapp * [rand([5000,1])*0.1];
end

g_k_ex = zeros(tsteps,1);
sra_pot_ex = zeros(tsteps,1);

if nargout > 2
    % Record REC (membrane voltage), Is (input currents), 
    % spk (spike raster), rs (firing rates) from all the units
    REC = zeros(N,tsteps);  % membrane voltage (in mV) values
%     Is = zeros(N, tsteps);  % input currents from the ext_stim
%     IPSCs = zeros(N, tsteps); % IPSC over time
%     spk = zeros(N, tsteps); % spikes
    rs = zeros(N, tsteps);  % firing rates
%     hs = zeros(N, tsteps); % filtered firing rates
end


I_noise = randn(size(Vm)) * noisescale;


for t=1:tsteps
    
    % calculate input currents from W and spks
    
    if mod(t,noiseduration) == 0
        I_noise = randn(size(Vm)) * noisescale;
    end
        
    
    if stp % NOTE STP is not set up right now
        % Short-term depression on the synaptic neurons, calculated in the
        % W_stp vector.
        
        RI = W  * (ri .* (1 - W_stp)) * syn_curr_multiplier;
        W_stp_idx(t) = W_stp(1);
    elseif ~isnan(Ksyn)
        RI = W * (Ksyn.*ri_slow * 0.9 + (1-Ksyn).*ri_fast) * syn_curr_multiplier;
    else
        % calculate input currents from W and spks
        RI = W * ri * syn_curr_multiplier;
    end
    
    if sra
        % Spike rate adaptation, include a potassium current to adapt
        % spiking rates to slow with more consecutive spiking.
        sra_pot = g_k .* (V_k - Vm);
        RI = RI + g_k .* (V_k - Vm);
        g_k_ex(t) = g_k(1);
        sra_pot_ex(t) = sra_pot(1);
    end
    
    if size(Iapp,2) > 1
        I = Iapp(:,t)*1e-2 + RI + I_noise;
    else
        I = Iapp*1e-2 + RI + I_noise; % Calculate synaptic currents
    end

    
    % Update membrane voltages using the update_LIF function
    [Vm, spks] = update_LIF(Vm,I);
    
    
    
   
    
    % Refractory period:
    % update_LIF updates Vm based on spks, so if we make sure to keep spks
    % under control here (i.e. make it so that a neuron cannot spike twice
    % in a row) then we should be good. Note this should be BEFORE spk_mat
    % gets updated.
    spks(ref ~= 0) = 0; % Re-set spikes to 0 if any neuron is still in refractory period.
    Vm(ref ~= 0) = -0.065; % Clamp refractory voltage to reset potential
    ref(ref ~= 0) = ref(ref ~=0) + 1; % Add a counter of 1 to refractory period
    ref(spks == 1) = 1; % Add 1 to any neurons that spiked this period.
    ref(ref > ref_period) = 0; % Re-set neurons who have passed their refractory period.
   
    spk_mat(:,t) = spks; % Update spike matrix
    
    % Synaptic Dynamics: Double Exponential Function
    
    if ~isnan(Ksyn) % If Ksyn is specified, we keep track of fast and slow synapses
        si_fast = si_fast + dt * (-si_fast/tau_rise_fast + (1/(tau_rise_fast * tau_fall_fast) * spks));
        ri_fast = ri_fast + dt * (-ri_fast/tau_fall_fast + si_fast);
        
        si_slow = si_slow + dt * (-si_slow/tau_rise_slow + (1/(tau_rise_slow * tau_fall_slow) * spks));
        ri_slow = ri_slow + dt * (-ri_slow/tau_fall_slow + si_slow);
    else
        si = si + dt * (-si/tau_rise + (1/(tau_rise * tau_fall) * spks));
        ri = ri + dt * (-ri/tau_fall + si);
    end
    
%     example_neuron(1,t) = Vm(example_neuron_idx);
%     example_neuron(2,t) = si(example_neuron_idx);
%     example_neuron(3,t) = ri(example_neuron_idx);
%     example_neuron(4,t) = I(example_neuron_idx);
    
    if stp
        % Here we adjust the synaptic efficacy of each neuron when it
        % fires, so that firing "uses up" 10% (or whatever stp_frac is) of
        % the synaptic resources available at the site. 
        if stp_oldmode
            W_stp = W_stp + dt * (-W_stp/stp_tau) + spks * stp_frac; % not self-similar decreases
        else
            W_stp = W_stp + dt * (-W_stp/stp_tau) + spks .* (1-W_stp) * stp_frac;
        % Old STP
        end
%         
    end
    
    if sra
        % Here we use a g_k conductance for spike rate adaptation, it
        % undergoes an exponential decay with immediate upwards increment
        % when a spike occurs. 
        g_k = g_k + dt * (-g_k/tau_k) + delta_g_k * spks;
    end
    
    if nargout > 1
        REC(:,t) = Vm;
        
        if ~isnan(Ksyn)
            rs(:,t) = (Ksyn.*ri_slow + (1-Ksyn).*ri_fast);
        else
            rs(:,t) = ri;
        end
    end
    
end

% if sra
%     figure
%     subplot(1,2,1)
%     plot(g_k_ex)
%     subplot(1,2,2)
%     plot(sra_pot_ex)
% end

% if plot
%     figure
%     addpath('../../Useful Functions')
%     plotSpikeRaster(logical(spk_mat),'PlotType','vertline2');
%     
%     figure
%     subplot(4,1,1)
%     plot(example_neuron(1,:)); title('Membrane Voltage of ex neuron');
%     subplot(4,1,2)
%     plot(example_neuron(2,:)); title('Internal current variable given for ex neuron');
%     subplot(4,1,3)
%     plot(example_neuron(3,:)); title('Outward current filter variable for ex neuron');
%     subplot(4,1,4)
%     plot(example_neuron(4,:)); title('Synaptic current applied to ex neuron');
%     
%     if stp
%         figure
%         plot(W_stp_idx)
%     end
%     if sra
%         figure
%         subplot(1,2,1)
%         plot(g_k_ex)
%         subplot(1,2,2)
%         plot(sra_pot_ex)
%     end
%     
% end

if ~isnan(p.Results.init_state),rng shuffle,end
    
end