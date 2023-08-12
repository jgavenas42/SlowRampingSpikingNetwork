function [spk_mat, rs, REC] = snap_singletrial_dimensionless(W,varargin)

% Simulate a single continuous trial of a LIF network

%% Parsing input

% Default input values
Iapp = 100; % Applied current in pA. Should be scalar, vector size Nx1, vector size (1xtsteps) or matrix size (Nxtsteps)
tsteps = 1200;
noisescale = 0; % standard deviation of membrane noise (resting -0.07)
noiseduration = 10; % how many timesteps during which to update the noise

curr_type = 'same'; % Type of current to administer. Can be 'same','supra', or 'sub'

tau_rise = .002;   % Synaptic rise time
tau_fall = .02;  % Synaptic fall time

% Alternative to specifying the tau_rise and tau_fall, specifying Ksyn will
% make it so that the function maintains activity of an average of fast and
% slow synapses, and Ksyn determines how they are averaged.
Ksyn = 0.0; 

p = inputParser;
addParameter(p,'Iapp',Iapp,@isnumeric);
addParameter(p,'tsteps',tsteps,@isnumeric);

addParameter(p,'noisescale',noisescale,@isnumeric);
addParameter(p,'noiseduration',noiseduration,@isnumeric);

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

if ~isnan(p.Results.init_state),rng(p.Results.init_state),end

Ksyn = p.Results.Ksyn;

curr_type = p.Results.curr_type;

N = size(W,1);
dt = 0.001;


% Preallocation of memory
spk_mat = zeros(size(W,1),tsteps);


Vm = zeros(size(W,1),1); % For keeping track of membrane potential
Vm = Vm + rand(size(Vm)); % randomly initialize membrane potential
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

%% If we want to give a slightly different size network
if strcmp(curr_type,'supra')
    Iapp = Iapp + Iapp * [rand([N,1])*0.1];
elseif strcmp(curr_type,'sub')
    Iapp = Iapp - Iapp * [rand([N,1])*0.1];
end



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
        
    
    if ~isnan(Ksyn)
        RI = W * (Ksyn.*ri_slow + (1-Ksyn).*ri_fast);
    else
        % calculate input currents from W and spks
        RI = W * ri;
    end
    
    
    if size(Iapp,2) > 1
        I = Iapp(:,t) + RI + I_noise;
    else
        I = Iapp + RI + I_noise; % Calculate synaptic currents
    end

    
    % Update membrane voltages using the update_LIF_dimensionless function
    [Vm, spks] = update_LIF_dimensionless(Vm,I);
    
   
    
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
    
    
    % If asked for, record synaptic activity AND 
    if nargout > 1
        REC(:,t) = Vm;
        
        if ~isnan(Ksyn)
            rs(:,t) = (Ksyn.*ri_slow + (1-Ksyn).*ri_fast);
        else
            rs(:,t) = ri;
        end
    end
    
end


if ~isnan(p.Results.init_state),rng shuffle,end
    
end