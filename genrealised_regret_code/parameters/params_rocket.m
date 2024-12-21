% Copyright (c) 2023, ETH Zurich, 
% Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, 
% Institute for Dynamic Systems and Control, D-MAVT
% All rights reserved.

function params = params_rocket(solver_SDP)
switch nargin
    case 0
        solver = '';
        
    case 1
        
    otherwise
        error('Wrong number of inputs!')
end

%% Simulation Parameters
params.sim.nrSteps = 25;
params.sim.nrTraj = 100;
params.sim.x_0 =  [0.7; 0.7;
                   0; 0;
                   0*pi/180; 0*pi/180];
               
%% System Parameters
% system dimensions
n = 6; m = 2;
params.sys.n = n;
params.sys.m = m;

% linear system
% states: y, z, v_y, v_z, theta, omega
% inputs: u_y, u_z
T_s = 0.1; %sampling time
mass = 1; %mass
l = 1; %length
J = 1; %moment of inertia
g = 9.81; %gravity

params.sys.A=[1, 0, T_s, 0, 0, 0; ...
              0, 1, 0, T_s, 0, 0; ...
              0, 0, 1, 0, 0, 0; ...
              0, 0, 0, 1, 0, 0; ...
              0, 0, 0, 0, 1, T_s; ...
              0, 0, 0, 0, 0, 1];
params.sys.B=[0, 0; ...
              0, 0; ...
              T_s/mass, 0; ...
              0, T_s/mass; ...
              0, 0; ...
              T_s*(J\l), 0];

% state constraints
params.sys.H_x = [eye(n);-eye(n)];
constraint_vals = [5;
                  5;
                  5;
                  5;
                  45*pi/180;
                  45*pi/180];
params.sys.h_x = [constraint_vals; constraint_vals]; 

% input constraints
u_max = 2*mass*g; % max pos. thrust deviation from eq. trust mass*g
u_min = 0.5*mass*g; % max neg. thrust deviation from eq. trust mass*g
delta_max = 30*pi/180; %max thrust angle

params.sys.H_u = [1,-tan(delta_max);-1,-tan(delta_max);0,1;0,-1];
params.sys.h_u = [mass*g*tan(delta_max); mass*g*tan(delta_max); u_max; u_min ]; 

% noise descriptions
params.sys.H_w = [eye(n);-eye(n)];
params.sys.noise_bounds = [3e-2;
                3e-2;
                3e-2;
                3e-2;
                7e-1*pi/180;
                7e-1*pi/180];
params.sys.h_w = [params.sys.noise_bounds; params.sys.noise_bounds];

%get all vertices
params.V = params.sys.noise_bounds'.*((-1).^ff2n(6));

%minimum volume ellipse inscribing the polytope
P = sdpvar(n,n);
objective = -logdet(P);
constraints = [P>=0];
for i=1:size(params.V,1)
    constraints = [constraints, params.V(i,:)*P*params.V(i,:)'<=1];
end
optimize(constraints, objective, sdpsettings('solver',solver_SDP,'verbose',0));
params.sys.P = value(P);

% max disturbance energy bound 
params.sys.omega = min(svd(params.sys.P))\params.sim.nrSteps;

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