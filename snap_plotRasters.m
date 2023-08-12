addpath(genpath('./utils'))

% Local
% fpath = './data_th/Ksyn_0.75/';

% Drive
fpath = '/Volumes/GDrive/snap/data_th/Ksyn_0/';

files = dir([fpath,'*.mat']);
% t_fid = fopen([fpath,'filesplotted2.txt'],'a+');
% t_A = fileread([fpath,'filesplotted2.txt']);

for i=1:length(files)
    t_fname = files(i).name;
    temp = split(files(i).name,'.mat');
    fname = temp{1};
    load([fpath,t_fname],'spks','N','tsteps')
    close(gcf)

    spks_trial = zeros(N,tsteps);

    for ii=1:N
        spks_trial(ii,:) = spks{ii}(1,:);
    end

    figure
    plotRaster(spks_trial);
    set(gca,'FontSize',16)
    xlabel('Time (ms)','FontSize',20)
    ylabel('Neuron #','FontSize',20)
    saveas(gcf,[fpath,fname,'/ExTrialRaster'],'epsc')
end
close all