function params = params_stab()
%% Simulation Parameters
params.sim.nrSteps = 30;
params.sim.x_0 =[3;1;3];
params.sim.nrTraj = 1000;

%% System Parameters
% system dimensions
n = 3; m = 2;
params.sys.n = n;
params.sys.m = m;

params.sys.rho=0.7;
% linear system
params.sys.A=params.sys.rho*[0.7, 0.2, 0; 0.3, 0.7, -0.1; 0, -0.2, 0.8];
params.sys.B=[1, 0.2; 2, 0.3; 1.5, 0.5];

% state constraints
params.sys.H_x = [1,0,0; -1,0,0; 0,1,0; 0,-1,0;0,0,1;0,0,-1];
params.sys.h_x = 3*ones(size(params.sys.H_x,1),1); 
% input constraints
params.sys.H_u = [1,0;-1,0;0,1;0,-1];
params.sys.h_u = 2*ones(size(params.sys.H_u,1),1); 
% noise description
params.sys.H_w = [1,0,0; -1,0,0; 0,1,0; 0,-1,0;0,0,1;0,0,-1];
params.sys.h_w = ones(size(params.sys.H_w,1),1);
params.sys.P = 1/n*eye(n);

% max disturbance bound
params.sys.omega = params.sys.n*params.sim.nrSteps;

%% Control Parameters
params.ctrl.T = params.sim.nrSteps;
params.ctrl.Sigma = eye(n);
params.ctrl.Q=eye(n);
params.ctrl.R=eye(m);
params.ctrl.x_0 = params.sim.x_0;
%Noise bound 
params.ctrl.P = params.sys.P;
%% Plot Parameters
params.plot.show = true;
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 1;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end