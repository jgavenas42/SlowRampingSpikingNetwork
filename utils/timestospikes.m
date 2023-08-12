function mat_spikes = timestospikes(fcn_spk)
    
% This function is designed to be used with SimLIFNet. When simulating a
% network of neurons, SimLIFNet outputs a cell of spike times rather than a
% matrix of zeros and ones corresponding to different time points. This
% function takes that cell of spike times and turns it into such a matrix.

% Note that the spike times are in seconds, so we multiply
% num_simtime by 1000 in order to make it milliseconds.


temp_spk_size = size(fcn_spk); % Grab how many neurons were simulated.
fcn_num_neurons = temp_spk_size(1); 
fcn_num_simtime = ceil(max([fcn_spk{:}])); % Grab maximum spike time.
mat_spikes = zeros(fcn_num_neurons,fcn_num_simtime*1000);

for i=1:fcn_num_neurons
    temp_vec = round(fcn_spk{i}.*1000);
    if isempty(temp_vec)
        continue
    end
    for j=1:length(temp_vec)
        mat_spikes(i,temp_vec(j))=1;
    end
end

        
        
        
        
        