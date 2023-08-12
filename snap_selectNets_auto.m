%% Main script for selecting connectivity matrices, etc. from selected data

t_fpath = '/Volumes/GDrive/snap/data_th_ex/noisescale_0.5/Ksyn_1/';

files = dir([t_fpath,'*.mat']);
t_fid = fopen([t_fpath,'filesselected_auto.txt'],'a+');
t_A = fileread([t_fpath,'filesselected_auto.txt']);

file_counter = 0;
base_period = 50; % how many continuous 50 ms bins "sub-threshold" (20 = 1 sec)
th_period = 1; % how many continuous 50 ms bins "above threshold"


for i=1:length(files)
    t_fname = files(i).name;
    t_idx = strfind(t_A,t_fname);
    
    searchstring = '07-27';
    
    if contains(t_fname,searchstring) && isempty(t_idx)
        
        disp(['Showing ',t_fname])
        load([t_fpath,t_fname],'ex_cluster_FRs','ex_tsteps','ex_Ntrials')
        
        [t_max,max_clust] = max(mean(ex_cluster_FRs,[1 3]));
        
        % Select only the cluster with maximum firing rate
        max_cluster_FRs = squeeze(ex_cluster_FRs(:,max_clust,:));
        
        % Normalize the firing rate similar to Fried et al
        max_Norm_FRs = (max_cluster_FRs - mean(max_cluster_FRs(:,1:10),2)) ./ (max(max_cluster_FRs,[],2));
        
        % Binarize to 0.5 normalized value, also similar to Fried et al
        % findings
        max_Norm_FRs_binary = max_Norm_FRs > 0.5;
        
        
        % What we want is a vector that has "baseline" activity for at
        % least 2 seconds and then hits some threshold... so create a
        % vector and match an appropriate vector to the binary vector.
        fits = zeros(ex_Ntrials,1);
        
        for k=1:ex_Ntrials
            fits(k) = ~isempty(strfind(max_Norm_FRs_binary(k,:),[repmat([0],base_period,1); repmat([1],th_period,1)]'));
        end
        
        if sum(fits) >= length(fits)/2
            fprintf(t_fid,[t_fname,'\n']);
            file_counter = file_counter + 1;
        end
        
        
    elseif ~contains(t_fname,searchstring)
        delete([t_fpath,t_fname]) 
    end
end

fclose(t_fid);