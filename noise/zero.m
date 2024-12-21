function [w] = zero(n,T)
%ZERO Creates zero noise
%   Outputs zero noise, i.e., no disturbance

%%% Parse inputs %%%
switch nargin
    case 2
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%

w = zeros(n,T);

end

