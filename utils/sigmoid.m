function [y] = sigmoid(x,varargin)

%% Sigmoid function for a vector x

% This function will take in a vector X, as well as optional parameters max,
% gain, and center, and run the usual logistic function on it.

% Inputs:
% x: a double-vector (required)
% gain: a scalar value that determines slope of the logistic function
% (default 1)
% center: a center value for the logistic (default 0)
% max: a maximum value for the logistic (default 1)

% Output:
% y: a double-vector output from logistic function


p = inputParser;
addParameter(p,'gain',1,@isnumeric);
addParameter(p,'center',0,@isnumeric);
addParameter(p,'max',1,@isnumeric);

parse(p,varargin{:})

y = (p.Results.max) ./ (1 + exp( - p.Results.gain * (x - p.Results.center)));

end






