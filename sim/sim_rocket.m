%% Setup

clc;
clear all;
close all;

% fix random number generator
rng(123);

% Select SDP Solver
solver_SDP = 'mosek';

%% get parameters and define system and controller

params = get_params('params_rocket',solver_SDP);


%% Synthesise controllers with state and input constraints
tic
ctrl_s{1} = NC_constr(params.sys, params.ctrl, solver_SDP);
toc
params.ctrl.O = ctrl_s{1}.O;

tic
ctrl_s{2} = H_2_constr(params.sys, params.ctrl, solver_SDP);
toc

tic
ctrl_s{3} = H_inf_constr(params.sys,params.ctrl, solver_SDP);
toc
tic
ctrl_s{4} = GenReg_constr(params.sys,params.ctrl,'dr',solver_SDP);
toc
tic
ctrl_s{5} = GenReg_constr(params.sys,params.ctrl,'cr',solver_SDP);
toc
tic
ctrl_s{6} = GenReg_PWB_constr(params.sys,params.ctrl,'dr',solver_SDP);
toc
tic
ctrl_s{7} = GenReg_PWB_constr(params.sys,params.ctrl,'cr',solver_SDP);
toc
%% simulate the closed-loop system for nrSteps time steps and nrTraj
% Define noise simulations
w_args{1} = {@gaussian_trunc_ell, {params.ctrl.T, zeros(params.sys.n,1), 5e-3*eye(params.sys.n) , params.sys.P},true};
w_args{2} = {@uniform_ell, {params.ctrl.T, params.sys.P},true};
w_args{3} = {@constant, {-params.sys.noise_bounds, params.sys.n, params.ctrl.T},false};
w_args{4} = {@sine, {params.sys.n, params.ctrl.T,-params.sys.noise_bounds, 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
w_args{5} = {@sawt, {params.sys.n, params.ctrl.T, -params.sys.noise_bounds, 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
w_args{6} = {@step, {params.sys.n, params.ctrl.T, params.sys.noise_bounds}, false};
w_args{7} = {@stair, {params.sys.n, params.ctrl.T, -params.sys.noise_bounds}, false};

%For all defined safe controllers and disturbance sims
for i=1:length(ctrl_s)
    for j=1:length(w_args)
        % If stochastic noise run nrTraj simulations
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
% All costs
for i=1:length(ctrl_s)-1
    for j=1:length(w_args)
        cost_comp_s(j,i)=cost_tot_s(end,i+1,j,1);
    end
end

%Normalize
for i=1:length(w_args)
    cost_comp_s(i,:)=cost_comp_s(i,:)/min(cost_comp_s(i,:));
end

% Controller names
for i=1:length(ctrl_s)-1
    ctrl_name_s{i}=ctrl_s{i+1}.name;
end

ctrl_name_s
cost_comp_s

%%
save(['RKT_T_' int2str(params.ctrl.T) '_x_' int2str(params.sim.x_0)' '_paper2'])