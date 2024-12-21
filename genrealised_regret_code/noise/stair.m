% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function [w] = stair(n,T,val)
%STAIR Creates stair noise
%   Outputs stair noise given dimension

%%% Parse inputs %%%
switch nargin
    case 3
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%

    w=val.*reshape([ones(n*(T-2*floor(T/3)),1); zeros(n*floor(T/3),1); ones(n*floor(T/3),1)],[n T]);
end

