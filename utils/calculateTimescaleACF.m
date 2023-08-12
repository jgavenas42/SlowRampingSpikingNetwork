function [tau_intrinsic, acf, inclusion] = calculateTimescale(varargin)

%% Calculate intrinsic timescale of a single neuron.
% Take a matrix of spiking activity across trials and calculate intrinsic
% timescale (timescale of autocorrelation decay) for a single neuron.

% Input:
% spk_mat:          binary spike matrix (n_trials X n_timepoints)

% Output:
% tau_intrinsic:    scalar-value time constant fitted on autocorrelation
% function. Will return NAN if doesn't pass inclusion criteria.
% acf:              aligned & averaged autocorrelation
% inclusion:        binary vector listing which inclusion criteria are
% violated.
% fitted:           the fitted model if we need more info.

% Outline: 
% 1) firingRate() to bin spiking activity into 50 ms bins
% 2) corr() to calculate cross-correlation between columns of firing rate
% matrix
% 3) alignMat() to take this diagonalized matrix and line up nonzero
% entries into left diagonal
% 4) Average aligned matrix
% 5) fit decaying autocorrelation function to obtain tau_intrinsic

% Inclusion Criteria:
% If a neuron does not pass one of the below four inclusion criteria,
% tau_intrinsic returns nan. Note that these checks might speed things up
% if fitting is slow, so check 1), then 4), then fit, then check 2) and 3).
% 1) neuron needs to have a minimum average firing rate of 1 Hz
% (experimental) or 2.5 Hz (simulated).
%       TO DO: implement experimental
% 2) tau needs to be between 0 and 500 ms. Removed this b/c some taus were
% higher but fine
% 3) A > 0 (amplitude of ACF)
% 4) first decrease in acf before 150 ms. 
% 5) the parameter search for tau_intrinsic converges/doesn't throw an
% error.


FR = nan;
spk_mat = nan;

p = inputParser;
addParameter(p,'FR',FR,@isnumeric);
addParameter(p,'spk_mat',spk_mat,@isnumeric);


parse(p,varargin{:})
FR = p.Results.FR;
spk_mat = p.Results.spk_mat;





%% Calculate firing rate

% firingRate takes a var (neurons, trials) X time binary spike matrix and
% returns a var X time/bin_size firing rate or spike count matrix.

if isnan(FR)
    if ndims(spk_mat) > 1
        spk_mat = squeeze(spk_mat);
    end
    % We use firing rate rather than spike count, because
    FR = firingRate(spk_mat','count',true,'bin_size',50);
    
end

inclusion = zeros(5,1);
if mean(FR,'all') < 2.5
    inclusion(1) = 1;
    tau_intrinsic = nan;
end

%% Calculate (auto)correlation
% As opposed to calculateTimescale, calculateTimescaleACF uses the standard
% ACF function to calculate autocorrelation

acf = autocorrFcn(FR',12,false);

acf_diff = diff(acf(1:end));

acf_diff_idx = find(acf_diff < 0); 

% Condition on first index where there is a decrease in the first 150
% milliseconds
if isempty(acf_diff_idx)
    inclusion(4) = 1;
    tau_intrinsic = nan;
    inclusion(2:3) = nan;
    acf = acf(1:12);
    return 
elseif (acf_diff_idx(1) > 3)
    inclusion(4) = 1;
    tau_intrinsic = nan;
    inclusion(2:3) = nan;
    acf = acf(1:12);
    return 
end

%% Fit autocorrelation function

% Create time lags vector for fitting using acf_diff_idx(1), i.e. the first
% decrease in autocorrelation.
timelags = acf_diff_idx(1)*50:50:600;

% Nonlinear least squares, use [1 1 1] as starting point because
% tau_intrinsic cannot equal zero (because we divide by it)
fo = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[0.5 1 50],...
               'Algorithm','Levenberg-Marquardt');%,...
%                'Lower',[-500 -500 -8000],...
%                'Upper',[500 500 8000]);
ft = fittype('A*(exp(-(x)/tau)+B)','options',fo);

% Fit timelags against acf from the first decrease until the 12th entry
% (corresponding to 600 milliseconds divided by 50.
try
    [fitted,gof,output] = fit(timelags',acf(acf_diff_idx(1):12),ft);
catch
    tau_intrinsic = nan;
    inclusion(5) = 1;
    acf = acf(1:12);
    return
end
    
tau_intrinsic = fitted.tau;

% Use lsqnonlin instead
% fun = @(x)x(1)*(exp(-timelags * x(2))+x(3))-acf(acf_diff_idx(1):12);
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
% x0 = [0.5 0.02 1];
% x = lsqnonlin(fun,x0,[],[],options);
% tau_intrinsic = x(2);

acf = acf(1:12);

if (tau_intrinsic < 0) || (tau_intrinsic > 500)
    tau_intrinsic = nan;
    inclusion(2) = 1;
end

if fitted.A < 0
% if x(1) < 0
    tau_intrinsic = nan;
    inclusion(3) = 1;
end


end