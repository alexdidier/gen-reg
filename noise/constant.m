function [w] = constant(val,n,T)
%CONSTANT Creates constant noise
%   Outputs constant noise given the value and dimension

%%% Parse inputs %%%
switch nargin
    case 1
        w=val;
    case 3
        w = val.*ones(n,T);
    otherwise
        error('Wrong number of inputs or val dimension!')
end
%%%%%%%%%%%%%%%%%%%

end

