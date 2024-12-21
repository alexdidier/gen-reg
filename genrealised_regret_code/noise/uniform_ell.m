% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function [w] = uniform_ell(T, P)
%UNIFORM Creates uniform noise in a given range
%   Samples a random disturbance vector uniformly from a given ellipsoid


%%% Parse inputs %%%
switch nargin
    case 2
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%

for i=1:T
    w(:,i) = chol(P,'L')\(2*rand(size(P,1),1)-1);
    while w(:,i)'*P*w(:,i)>1
        w(:,i) = chol(P,'L')\(2*rand(size(P,1),1)-1);
    end
end
end

