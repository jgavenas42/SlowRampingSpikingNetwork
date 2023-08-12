%% Main script for selecting connectivity matrices, etc. from selected data

t_fpath = '/Volumes/GDrive/snap/data_th/noisescale_0.5/Ksyn_0.75/';

files = dir([t_fpath,'*.mat']);
t_fid = fopen([t_fpath,'filesselected.txt'],'a+');
t_A = fileread([t_fpath,'filesselected.txt']);

for i=1:length(files)
    t_fname = files(i).name;
    t_idx = strfind(t_A,t_fname);
    
    searchstring = '07-27';
    
    if contains(t_fname,searchstring) && isempty(t_idx)
        
        disp(['Showing ',t_fname])
        load([t_fpath,t_fname],'ex_cluster_FRs','ex_tsteps','ex_Ntrials')
        
        f=figure('units','normalized','outerposition',[0 0 1 1]);
        for j=1:ex_Ntrials
            subplot(ex_Ntrials,1,j)
            plot((1:50:ex_tsteps)/(1000),squeeze(ex_cluster_FRs(j,:,:))'), hold on
            plot((1:50:ex_tsteps)/(1000),mean(squeeze(ex_cluster_FRs(j,:,:)),1),'k','LineWidth',2)
            xlabel('time (s)'),ylabel('Firing Rate (Hz)')
            title('FRs on each cluster + avg')

%             subplot(i,2,2)
%             plot((1:ex_tsteps)/1000,movmean(cluster_PSPs',50)), hold on
%             plot((1:ex_tsteps)/1000,movmean(mean(cluster_PSPs,1),50),'k','LineWidth',2)
%             xlabel('time (s)'),ylabel('mV (post synaptic potentials)'),title('PSPs on each cluster + avg')
        end
        x = input('Does it look good? '); % Will ask if the input is good, if so push 1, if not push 0.
        
        if x == 1
            fprintf(t_fid,[t_fname,'\n']);
        end
        
        close(f)
        
    elseif ~contains(t_fname,searchstring)
        delete([t_fpath,t_fname]) 
    end
end

fclose(t_fid);