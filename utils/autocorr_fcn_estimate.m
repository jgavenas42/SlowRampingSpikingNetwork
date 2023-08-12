function [m epoch xepoch] = autocorr_fcn_estimate(singleTrace,epochLen,nSamplings)
%
% function [m epoch xepoch] = autocorr_fcn_estimate(singleTrace,epochLen,nSamplings)
% 
% This will be to estimate the autocorr fcn for a single trace of signal data
%

epoch = zeros(nSamplings,epochLen);
xepoch = zeros(nSamplings,epochLen*2-1);

startPoints = ceil(rand(nSamplings,1)*(length(singleTrace)-epochLen));
for i=1:nSamplings
    epoch(i,:) = singleTrace(startPoints(i):startPoints(i)+(epochLen-1));
    xepoch(i,:) = xcorr(epoch(i,:),'coeff');
end

% this is the left-hand side of the xcorr fcn
m = mean(xepoch(:,1:epochLen))';

