
addpath(genpath('./utils'))

fpath = '/Users/gavenas/Desktop/Research/snap/fw_fix copy/';

fix = 'yes';

foldernames = dir([fpath,'sub*']);

savepath = ('data_real/');

correct_trials = xlsread([fpath,'correct_trials.xls']);

for i=1:length(foldernames)
    
    sinfo = foldernames(i).name;
    
    filenames = dir([fpath,sinfo,'/*.mat']);
    
    
    
    for j=1:length(filenames)
        
        finfo = filenames(j).name;
        
        load([fpath,sinfo,'/',finfo])
        
        data = [];
        
        data.BaselineAlignedSpikes = BaselineAlignedSpikes;
        data.StopAlignedSpikes = StopAlignedSpikes;
        data.WAlignedSpikes = WAlignedSpikes;
        
        
        t_idx = strfind(sinfo,'s');
        data.subj = str2double(sinfo(5:t_idx(2)-1));
        data.session = str2double(sinfo(end));
        
        
            
        
        t_info = split(finfo,'_');
        
        data.channel = str2double(t_info{1}(8:end));
        data.cluster = str2double(t_info{2}(8:end));
        data.region = t_info{3};
        
        data.ntrials = size(BaselineAlignedSpikes,1);
        
        % Do some fixing to see if the file has the right number of trials.
        % This 'correct trials' sheet was created manually by looking at
        % the data, which had been mixed, and figuring out which subjects &
        % sessions had certain number of trials.
        if strcmp(fix,'yes')
            ix_row = find(correct_trials(:,1)==data.subj & correct_trials(:,2)==data.session);
            
            if correct_trials(ix_row,3) ~= data.ntrials
                delete([fpath,sinfo,'/',finfo])
                continue
            end
        end
        
        data.time = linspace(-3.5,0,size(StopAlignedSpikes,2));
        
        idxs = ((data.time >= -3) & (data.time < -2));
        baseFRs = mean(StopAlignedSpikes(:,idxs),2)*1000;
%         baseFRs = mean(BaselineAlignedSpikes,2)*1000;
        data.baseFRs = baseFRs;
        
%         stopFRs = mean(WAlignedSpikes(:,3101:end),2)*1000;
        stopFRs = mean(StopAlignedSpikes(:,3101:end),2)*1000;

        data.stopFRs = stopFRs;
        
        data.FRs = firingRate(StopAlignedSpikes);
%         data.FRs = firingRate(WAlignedSpikes);

        
        % Rank-sum test to classify
        [P,H,STATS] = ranksum(stopFRs,baseFRs,'alpha',0.01);
        
        data.p = P;
        data.zval = STATS.zval;
        if isnan(H) || nanmean(data.FRs,'all') < 0.5
            data.type = 'bad';
        elseif ~H
            data.type = 'nochange';
        elseif (STATS.zval > 0)
            data.type = 'inc';
        elseif (STATS.zval < 0)
            data.type = 'dec';
        end
        
        % Sigmoid fitting
        % Fitting Sigmoid
        fB = num2str(mean(baseFRs));
        fW = num2str(mean(stopFRs));
        
        A = [];
        for bi=1:data.ntrials
            A = [A diff(find(StopAlignedSpikes(bi,:)==1))];
        end
        data.BI=mean(A<20);

        if strcmp(data.type,'inc')
            model = strcat('(',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
            
            fo = fitoptions('Method','NonlinearLeastSquares',...
                           'Algorithm','Levenberg-Marquardt');
            ft = fittype(model,'options',fo);

            frs = nanmean(data.FRs,1);
            time_fr = linspace(-3.5,0,length(frs));
            idxs = find((time_fr > -2.5));
            [fitted,gof,output] = fit(time_fr(idxs)',smoothdata(frs(idxs),'gaussian',4)',ft);
            data.sigmoid = fitted;
            data.sigmoid_gof = gof;
            data.gain = fitted.gain;

            plot(fitted,time_fr(idxs),smoothdata(frs(idxs),'gaussian',4))
            xlim([-2.5 0])
            legend
        elseif strcmp(data.type,'dec')
            model = strcat('(-',fW,'/(1 + exp(-(x - t0)/gain)))+',fB);
            fo = fitoptions('Method','NonlinearLeastSquares',...
                           'Algorithm','Levenberg-Marquardt');
            ft = fittype(model,'options',fo);

            frs = nanmean(data.FRs,1);
            time_fr = linspace(-3.5,0,length(frs));
            idxs = find((time_fr > -2.5));
            [fitted,gof,output] = fit(time_fr(idxs)',smoothdata(frs(idxs),'gaussian',4)',ft);
            data.sigmoid = fitted;
            data.sigmoid_gof = gof;
            data.gain = fitted.gain;

            plot(fitted,time_fr(idxs),smoothdata(frs(idxs),'gaussian',4))
            xlim([-2.5 0])
            legend
            
        end

        % Norm / heterogeneity analysis
        
        spks_sm = smoothdata(StopAlignedSpikes(:,501:end),2,'gaussian',400)*1000;
        spks_sm = spks_sm / max(spks_sm,[],'all');
        spks_sm_mean = nanmean(spks_sm,1);

        nCk = nchoosek(1:size(spks_sm,1),2);

        norms_this = 0;
        for kk=1:length(nCk)
        %         norms_this = norms_this + norm(spks_sm(kk)-spks_sm_mean);
            norms_this = norms_this + norm(spks_sm(nCk(kk,1))-spks_sm(nCk(kk,2)));

        end
        data.avgdist = norms_this / length(nCk);
        
        t_normFRs = (data.FRs - baseFRs);
        
        data.normFRs = t_normFRs ./ max(abs(t_normFRs'))';
        
        data.FRs = firingRate(StopAlignedSpikes);
        data.spks = StopAlignedSpikes;
        
        data.tau = calculateTimescale('spk_mat',data.BaselineAlignedSpikes);
        
        
        
        save([savepath,'subj',num2str(data.subj),'_ses',num2str(data.session),...
            '_chan',num2str(data.channel),'_clust',num2str(data.cluster),...
            '.mat'],'-struct','data')
        
    end
end

        
        
        
        
        
        
        