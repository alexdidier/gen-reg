function [w] = gaussian_trunc_poly(T, mean, covariance, H_w, h_w)
%GAUSSIAN Creates Gaussian distubance
%   Computes Gaussian disturbance based on noise mean and covariance

%%% Parse inputs %%%
switch nargin
    case 5
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%
for i=1:T
    w(:,i) = mean + chol(covariance,'L')*randn(size(mean));
    while any(H_w*w(:,i)>h_w)
        w(:,i) = mean + chol(covariance,'L')*randn(size(mean));
    end
end
end

