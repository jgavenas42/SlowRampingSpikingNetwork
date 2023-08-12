function firing_rates = firingRate(mat_spikes,varargin)

% Calculate firing rates from a spike matrix using nonoverlapping bins.

% Only ones and zeros should comprise the spike matrix, not spike times.

% Returns a matrix which is (size of matrix)/(bin size) containing the
% neural firing rates at each bin time.

% Resulting matrix: rows = neurons, columns = bins/time.

default_binsize = 50; % default 50 ms bins for 1 ms spike matrices.
default_units = "seconds"; % default unit seconds, used for calculating firing rates in Hz.
p = inputParser;
addParameter(p,'bin_size',default_binsize,@isnumeric);
addParameter(p,'units',default_units,@isstring);
addParameter(p,'count',false);

parse(p,varargin{:})

% Raise a warning if matrix and binsize are going to be cut off. Will round
% number of columns/bin size down and give that many firing rate values.
if(fix(size(mat_spikes,2)/p.Results.bin_size)~=size(mat_spikes,2)/p.Results.bin_size)
    warning("Matrix and bin size not compatible, spikes at end may be lost");
end    

% Create our empty firing rates vector, which we shall fill below.
firing_rates = zeros(size(mat_spikes,1),fix(size(mat_spikes,2)/p.Results.bin_size));

% Iterate through rows and then columns of firing_rates vector
for i=1:size(firing_rates,1)
    for j=1:size(firing_rates,2)
        
        % Counts the number of non-zero entries in each bin in mat_spikes
        firing_rates(i,j) = sum(mat_spikes(i,(1+p.Results.bin_size*(j-1)):(p.Results.bin_size*j)));
    end
end
if(p.Results.units ~= "seconds")
    warning("Firing rate units unclear");
elseif(p.Results.units=="seconds")
    if(not(p.Results.count))
        firing_rates=firing_rates*(1000/p.Results.bin_size);
    end
end

end


    
    
