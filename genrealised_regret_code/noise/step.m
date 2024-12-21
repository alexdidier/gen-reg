% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function [w] = step(n,T,val)
%STEP Creates step function noise
%   Outputs step function noise given and dimension

%%% Parse inputs %%%
switch nargin
    case 3
        
    otherwise
        error('Wrong number of inputs or val dimension!')
end
%%%%%%%%%%%%%%%%%%%
    w=val.*reshape([zeros(n*(T-floor(T/2)),1); ones(n*floor(T/2),1)],[n T]);
end

