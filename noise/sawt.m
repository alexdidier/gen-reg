function [w] = sawt(n,T,a,f,p)
%SAWT Creates sawtooth noise
%   Outputs sawtooth noise, with given frequency and phase

%%% Parse inputs %%%
switch nargin
    case 5
        
    otherwise
        error('Wrong number of inputs!')
end
%%%%%%%%%%%%%%%%%%%
t = 1:T;
for i=1:n
    w(i,:) = a(i)*sawtooth(2*pi*f(i)*(t-p(i)));
end
end

