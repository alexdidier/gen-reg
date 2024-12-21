function [x, u, cost, cost_tot,w] = sim_LQ(ctrl,w_f,w_arg,A,B,Q,R,x_0,T,time_avg)
%SIM_LQ Simulate system
%  Simulate a system given system matrices A, B, controller matrix K, cost
%  matrices Q, R, disturbance w, time steps T

switch nargin
    case 9
        time_avg = false;
        
    case 10
        
    otherwise
        error('Wrong number of inputs!');
        
end

n = size(A,2);
m = size(B,2);

x = zeros(n,T+1);
u = zeros(m,T+1);
cost = zeros(1,T+1);

x(:,1) = x_0;

w=w_f(w_arg{:});


for k=1:T
    if isprop(ctrl, 'nc')
        u(:,k) = ctrl.inp([x_0, w],k);
    else
        u(:,k) = ctrl.inp(x,k);
    end
    x(:,k+1) = A*x(:,k) + B*u(:,k) + w(:,k);
    cost(1,k) = x(:,k)'*Q*x(:,k)+u(:,k)'*R*u(:,k);
end

cost(1,T+1) = x(:,T+1)'*Q*x(:,T+1);

cost_tot=zeros(1,T+1);
cost_tot(1) = cost(1);

for i=2:T+1
    cost_tot(i) = cost_tot(i-1) +cost(i);
end

if time_avg
    for i=1:T+1
        cost_tot(i)=cost_tot(i)/i;
    end
end
end

