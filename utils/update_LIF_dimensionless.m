function [Vm, spks] = update_LIF_dimensionless(Vm,I)

% Function for updating membrane potential for LIF neurons
% Dimensionless version of update_LIF
%       I.E. resets to 0, fires at 1

% Inputs:

% Required:
% Vm:       Membrane voltage at time t-1
% adjacency matrix)
% Iapp:     Applied current (either scalar or vector of size Vm)

% Optional (from varargin): TO-DO
% Isyn_update:   Function for updating synaptic currents (e.g. conductance
% based decaying current, rise+fall current, etc.
% Isyn_params:   Parameters for updating synaptic currents.
% refractory:   Function for clamping voltages at reset potential for some
% time.
% refrac_params:    Parameters for refractory period.

% Outputs:
% Vm:       Updated membrane voltage at time t
% spks:     Which neurons spiked at time t



% if narging > 3: TO-DO
    

%% LIF Parameters
tau = 0.010;                % membrane time constant
E_L = 0;               % leak potential (also resting potential)
Vth = 1;               % threshold potential (to produce spike)
Vreset = 0;            % reset potential (post-spike)
% Cm = 100e-12;               % total membrane capacitance
% G_L = Cm/tau;               % total membrane conductance (leak conductance)
dt = 0.001;                 % Time step in milliseconds

%% Housekeeping

% Update the input current I with spiking activity, applied current, etc.
% CHANGE: handle the changes of the current I in the outer function and
% just pass the I for each trial here. 

% I = Iapp + W * spks;

% Reset spks only after calculating I
spks = zeros(size(Vm));

%% Update LIF


Vm = Vm + dt*(I + (E_L-Vm)/tau); %+ randn(size(Vm))*0.0001; % Dynamic Updating
spks(Vm > Vth) = 1;                 % Update spikes
Vm(Vm > Vth) = Vreset;              % Reset potential based on spiking

end

