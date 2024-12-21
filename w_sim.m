%% Setup

clc;
clear all;
close all;

for counter = 1:1
% fix random number generator
rng(123);

% Select SDP Solver
solver_SDP = 'mosek';
save_inter = false;
%% get parameters and define system and controller

load('RKT_T_25_x_110000_final3')
%params.sys.noise_bounds=params.sys.h_w(1:params.sys.n);
%params.sys.W = Polyhedron(params.sys.H_w,params.sys.h_w);
%w_args={};
params.sim.nrTraj = 100;
%params.sim.x_0=zeros(params.sys.n,1);
% params.ctrl.x_0=params.sim.x_0;
% params.sys.x_0=params.sim.x_0;
%% simulate the closed-loop system for nrSteps time steps and nrTraj
% Define noise simulations
w_args{1} = {@gaussian_trunc_ell, {params.ctrl.T, zeros(params.sys.n,1), 5e-3*eye(params.sys.n), params.sys.P}, true}; %inv(params.sys.P) , params.sys.P},true};
w_args{2} = {@uniform_ell, {params.ctrl.T, params.sys.P},true};
w_args{3} = {@constant, {-params.sys.noise_bounds, params.sys.n, params.ctrl.T},false};
w_args{4} = {@sine, {params.sys.n, params.ctrl.T,-params.sys.noise_bounds, 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
w_args{5} = {@sawt, {params.sys.n, params.ctrl.T, -params.sys.noise_bounds, 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
w_args{6} = {@step, {params.sys.n, params.ctrl.T, params.sys.noise_bounds}, false};
w_args{7} = {@stair, {params.sys.n, params.ctrl.T, -params.sys.noise_bounds}, false};

%%For all defined safe controllers and disturbance sims
for i=1:length(ctrl_s)
    for j=1:length(w_args)
         if w_args{j}{3}==true
            %Ensure same disturbances for every different controller
            rng(123)
            for k=1:params.sim.nrTraj
                [x_s(:,:,i,j,k), u_s(:,:,i,j,k), cost_s(:,i,j,k), cost_tot_s(:,i,j,k+1),w_s(:,:,i,j,k)] = sim_LQ(ctrl_s{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
            end
            cost_tot_s(end,i,j,1)=sum(cost_tot_s(end,i,j,2:end))/params.sim.nrTraj;
        else 
            [x_s(:,:,i,j,1), u_s(:,:,i,j,1), cost_s(:,i,j,1), cost_tot_s(:,i,j,1),w_s(:,:,i,j,1)] = sim_LQ(ctrl_s{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
        end
    end
end

%% Safe controllers
for i=1:length(ctrl_s)-1
    for j=1:length(w_args)
        cost_comp_s(j,i)=cost_tot_s(end,i+1,j,1);
    end
end

%Normalize
for i=1:length(w_args)
    cost_comp_s(i,:)=cost_comp_s(i,:)/min(cost_comp_s(i,:));
end

for i=1:length(ctrl_s)-1
    ctrl_name_s{i}=ctrl_s{i+1}.name;
end

ctrl_name_s
cost_comp_s

if cost_comp_s(1,1)==1
    break
end
end

%%
%save(['T_' int2str(params.ctrl.T) '_RKT_x_' int2str(params.sim.x_0)'])
