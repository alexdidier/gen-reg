function params = params_GenReg()
%% System Parameters
n = 4; m = 3;
params.sys.n = n;
params.sys.m = m;
params.sys.A=randn(n,n);
params.sys.B=randn(n,m);

params.sys.C=eye(n);
params.sys.D=zeros(n,m);

% state constraints
%params.sys.H_x = [1,0; -1,0; 0,1; 0,-1];
%params.sys.h_x = 10*[pi/4; pi/4; pi/3; pi/3]; %[rad]
% input constraints
%params.sys.H_u = [1;-1];
%params.sys.h_u = 10*[5;5]; %[1/s^2]
% noise description
%params.sys.H_w = [];
%params.sys.h_w = [];
params.sys.P = [5,4,3,2;4,6,2,1;3,2,7,2;2,1,2,5];

% max disturbance bound
params.sys.omega = 10;

%% Simulation Parameters
params.sim.nrSteps = 10;
params.sim.nrTraj = 1;
%params.sim.x_0 =[1;10];
params.sim.x_0=3*rand(n,1);
%% Control Parameters
params.ctrl.T = params.sim.nrSteps;
params.ctrl.Sigma = eye(n);
params.ctrl.Q=0.1*[5,4,3,2;4,6,2,1;3,2,7,2;2,1,2,5];%
params.ctrl.R=[5,4,3;4,6,2;3,2,7];%
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