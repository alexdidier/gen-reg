function [w] = sine(n,T,a,f,p)
%SINE Creates sinusoidal noise
%   Outputs sinusoidal noise, with given amplitude, frequency and phase

%%% Parse inputs %%%
switch nargin
    case 5
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%

for i=1:n
    for j=1:T
        w(i,j) = a(i)*sin(f(i)*j+p(i));
    end
end

end

