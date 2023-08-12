function p = plot_signal_ci(x,y,varargin)

p = inputParser;
p.addRequired('x');
p.addRequired('y')
p.addParameter('FigHandle',gcf,@isinteger);
p.addParameter('Color','black'); % specify nan for color & colors for exc & inh if separate 
% p.addParameter('N',400,@isinteger);
p.addParameter('LineWidth',2.5);
p.addParameter('FaceAlpha',0.2);

p.parse(x,y,varargin{:});


% This function will plot a mean signal with standard deviation from a set
% of 

N = size(y,1);                                      % Number of ?Experiments? In Data Set
yMean = nanmean(y);                                    % Mean Of All Experiments At Each Value Of ?x?
ySEM = nanstd(y)/sqrt(N);                              % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));   % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?

p2 = patch([x fliplr(x)], [yCI95(1,:)+yMean fliplr(yCI95(2,:)+yMean)],...
    p.Results.Color,'FaceAlpha',p.Results.FaceAlpha);
p2.EdgeColor= 'none';

hold on
p = plot(x, yMean,'Color',p.Results.Color,'LineWidth',p.Results.LineWidth);                                      % Plot Mean Of All Experiments



% hold off
% grid