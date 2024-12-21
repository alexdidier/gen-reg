function [w] = gaussian(T, mean, covariance)
%GAUSSIAN Creates Gaussian distubance
%   Computes Gaussian disturbance based on noise mean and covariance

%%% Parse inputs %%%
switch nargin
    case 3
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%
for i=1:T
    w(:,i) = mean + chol(covariance,'L')*randn(size(mean));
end
end

