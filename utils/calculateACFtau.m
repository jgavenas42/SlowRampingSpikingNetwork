function [tau, acf, lags] = calculateACFtau(spks,doplot)
%CALCULATEACFTAU Summary of this function goes here
%   Detailed explanation goes here

binsize = 50;
maxlag = 1000;

% lags = 1:1000;
lags = binsize:binsize:maxlag;
ix2use = lags < 1500;

acfs = nan(size(spks,1),size(lags(ix2use),2));

neuron2use = 1;
changeneuron = 1;


for ii=1:size(spks,1)
    
    
    x = spks(ii,:);
%     x = smoothdata(x,'gaussian',50)*1000;

    xx = reshape(x,binsize,[]);
    xx = sum(xx);
    
    if sum(x)*1000/length(x) < 2.5
        continue
    elseif changeneuron == 1;
        neuron2use = ii;
        changeneuron = 0;
    end
%     if mean(x) < 0.5,continue,end

%     YY = autocorrFcn(xx',max(lags),false);
    YY = autocorrFcn(xx',maxlag/binsize,false);

        
    acfs(ii,:) = YY(ix2use);
    
%     acfs
    
end

x = spks(neuron2use,:);

acf = nanmean(acfs,1);
sem = nanstd(acfs,1)/sqrt(sum(~isnan(acfs(:,1))));

lags = lags(ix2use);

acf_diff = diff(acf);

acf_diff_ix = find(acf_diff < 0);

if acf_diff_ix(1) > round(150 / binsize)
    tau = nan;
else
    try
        g = fittype('a+b*exp(-x/tau)');
        f = fit(lags(acf_diff_ix(1:end))',acf(acf_diff_ix(1:end))',g,'StartPoint',[0,1,100]);
        if f.tau > 1000 || f.tau < 0 || f.b < 0
            tau = nan;
        else
            tau = f.tau;
        end
    catch
        tau = nan;
    end
%         f = fit(lags(ix2use)',YY(ix2use),'exp1','StartPoint',[1,-0.1]);
end

if doplot==1 
        
    figure
    tileplot = tiledlayout(1,2);
    nexttile
    plot(x,'linewidth',3)
    xlabel('Time (ms)')
    ylabel('Example Neuron Spiking')
    nexttile
    plot(lags,acf,'k.','LineWidth',5)
    if ~isnan(tau)
        hold on
        xx = 0:1000;
        plot(xx,f.a + f.b * exp(-xx/f.tau),'LineWidth',3)
    end
    xlabel('Time Lag (ms)')
    ylabel('Autocorrelation & Fitted Exponential')
    title(sprintf("Autocorr Timescale = %0.2f",tau))

end

