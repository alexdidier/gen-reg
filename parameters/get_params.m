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