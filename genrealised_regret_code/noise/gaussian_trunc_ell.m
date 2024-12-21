% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function [w] = gaussian_trunc_poly(T, mean, covariance, P)
%GAUSSIAN Creates Gaussian distubance
%   Computes Gaussian disturbance based on noise mean and covariance

%%% Parse inputs %%%
switch nargin
    case 4
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%
for i=1:T
    w(:,i) = mean + chol(covariance,'L')*randn(size(mean));
    while any(w(:,i)'*P*w(:,i)>1)
        w(:,i) = mean + chol(covariance,'L')*randn(size(mean));
    end
end
end

