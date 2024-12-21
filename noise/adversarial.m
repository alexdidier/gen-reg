function [w] = adversarial(Phi, C, n, T, x_0, constr, P)
%ADVERSARIAL Creates adversarial noise given system response
%   Outputs adversarial noise given a system response and noise constraints

%%% Parse inputs %%%
switch nargin
    case 7

    otherwise
        error('Wrong number of inputs or val dimension!')
end
%%%%%%%%%%%%%%%%%%%

w=sdpvar(n*T,1);
objective = -[x_0;w]'*(Phi'*C*Phi)*[x_0;w];

switch constr
    case 'PWB'
        constraints = [];
        for i=1:T
            P_i=zeros(n*T);
            P_i((i-1)*n+1:i*n,(i-1)*n+1:n*T)=P;
            constraints=[constraints, w'*P_i*w<=1];
        end
        
    case 'BE'
        constraints=[w'*P*w<=1];
        
    case 'Poly'
        constraints = [kron(eye(T),P)*w<=ones(T*size(P,1),1)];
    
    otherwise
        error('Not a defined constraint type!');
end
disp('Solving...')
tic
optimize(constraints, objective, sdpsettings('verbose',1,'solver','fmincon'));
toc
w = reshape(value(w),[n T]);
end

