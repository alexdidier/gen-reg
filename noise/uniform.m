function [w] = uniform(T, W)
%UNIFORM Creates uniform noise in a given range
%   Samples a random disturbance vector uniformly from a given polytope
%   using MPT3

%%% Parse inputs %%%
switch nargin
    case 2
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%

for i=1:T
    w(:,i)=W.randomPoint();
end
end

