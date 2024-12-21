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

