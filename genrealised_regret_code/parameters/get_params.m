% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function params = get_params(arg,solver_SDP)
%%%%%%%%%%%%%%
switch nargin
    case 1
        params = eval(arg);
        
    case 2
        params = eval([arg,'(solver_SDP)']);
        
    otherwise
        error('Wrong number of inputs!')
end

%%%%%%%%%%%%%%
end