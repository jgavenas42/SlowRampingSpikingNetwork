function plotRaster(spikes,varargin)

%% Plot a raster given a binary spike matrix

% Inputs:
% Spks - binary spike matrix (must be 2D)
% varargin - used for plotting


%%
p = inputParser;
p.addRequired('spikes');
p.addParameter('FigHandle',gcf,@isinteger);
p.addParameter('Color','black'); % specify nan for color & colors for exc & inh if separate 
p.addParameter('ColorSep',nan);
p.addParameter('Color_exc','red');
p.addParameter('Color_inh','blue');
% p.addParameter('N',400,@isinteger);
p.addParameter('perc_exc',0.8,@isinteger);
p.addParameter('MarkerSize',0.5)

p.parse(spikes,varargin{:});


[A,B] = find(spikes == 1);

perc_exc = p.Results.perc_exc;
N = size(spikes,1);
MS = p.Results.MarkerSize;

if ~isnan(p.Results.ColorSep)
    plot(B(A <= perc_exc * N),A(A<=perc_exc*N),'.','Color',p.Results.Color_exc,'MarkerSize',MS)
    hold on
    plot(B(A > perc_exc * N),A(A>perc_exc*N),'.','Color',p.Results.Color_inh,'MarkerSize',MS)
else
    plot(B,A,'.','Color',p.Results.Color,'MarkerSize',MS)
end
ylim([0,size(spikes(~isnan(spikes(:,1,1)),:,:),1)+1])


