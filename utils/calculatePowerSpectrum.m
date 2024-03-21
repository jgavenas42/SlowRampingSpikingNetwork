function [fexp,Pxx,frex] = calculatePowerSpectrum(spks,doplot)
%CALCULATEPOWERSPECTRUM Summary of this function goes here
%   Detailed explanation goes here


if nargout > 1
    fexps = nan(size(spks,1),1);
end

Fs = 1000; % 1000 Hz sampling rate, shouldn't change...
T = 1/Fs;             % Sampling period       
L = size(spks,2);             % Length of signal, shouldn't change
t = (0:L-1)*T;        % Time vector

frex = Fs/L*(0:L-1);
ix2use = frex > 0 & frex < 20;

Pxxs = nan(size(spks,1),size(frex(ix2use),2));


for ii=1:size(spks,1)
    
    
    x = spks(ii,:);
    x = smoothdata(x,'gaussian',400)*1000;
    
    if mean(x) < 0.5
        continue
    end
    
    Y = fft(x);
    YY = abs(Y);
    
    
    Pxxs(ii,:) = YY(ix2use);
    


    if nargout > 1
        
        p = polyfit(log(frex(ix2use)),log(YY(ix2use)),1);
        
        fexps(ii) = p(1);
        
    end
    
    if doplot==1 
        
        doplot = 0;
        figure
        tileplot = tiledlayout(3,1);
        nexttile
        plot(t,x,'linewidth',3)
        xlabel('Time (ms)')
        ylabel('Spiking')
        nexttile
        YY = abs(Y);
        plot(frex(ix2use),YY(ix2use),'linewidth',3)
        xlabel('Frequency (Hz)')
        ylabel('|fft(spiking)|')
        nexttile
        xx = log(frex(ix2use));
        plot(xx,log(YY(ix2use)))
        
        
        xlabel('log(Frequency (Hz))')
        ylabel('log(|fft(spiking)|)')
        
        if nargout > 1
            hold on
            plot(xx,p(1)*xx+p(2),'r','linewidth',2)
        end
    end
end

if nargout > 1
    fexp = nanmean(fexps);
end

Pxx = nanmean(Pxxs,1);

frex = frex(ix2use);


    

