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

