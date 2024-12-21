%% Setup

%clc;
clear all;
close all;

% fix random number generator

% Select SDP Solver
solver_SDP = 'mosek';

%% get parameters and define system and controller

params = get_params('params_Rand');
params.sys.A
%ctrl{1} = NC(params.sys, params.ctrl);
%ctrl{2} = H_2(params.sys, params.ctrl, solver_SDP);
%ctrl{3} = H_inf(params.sys,params.ctrl, solver_SDP);

%ctrl{4} = GenReg(params.sys, params.ctrl,'dr', solver_SDP);
%ctrl{5} = GenReg(params.sys, params.ctrl, 'cr', solver_SDP);

ctrl{1} = GenReg_PWB(params.sys,params.ctrl,'dr',solver_SDP);
ctrl{2} = GenReg_PWB(params.sys,params.ctrl,'cr',solver_SDP);

ctrl{3} = GenReg_FBS(params.sys,params.ctrl,'dr',solver_SDP);
ctrl{4} = GenReg_FBS(params.sys,params.ctrl,'cr',solver_SDP);

%ctrl{10} = GenReg_PWB_constr(params.sys,params.ctrl,'dr',solver_SDP);
%ctrl{11} = GenReg_PWB_constr(params.sys,params.ctrl,'cr',solver_SDP);

%ctrl{12} = GenReg_FBS_constr(params.sys,params.ctrl,'dr',solver_SDP);
%ctrl{13} = GenReg_FBS_constr(params.sys,params.ctrl,'cr',solver_SDP);

% %%
% a=rand(params.sys.n*(params.ctrl.T+1));
% ctrl_aPWB = GenReg_PWB(params.sys,params.ctrl,a*a'+0.1*eye(params.sys.n*(params.ctrl.T+1)),solver_SDP);
% ctrl_aFBS = GenReg_FBS(params.sys,params.ctrl,a*a'+0.1*eye(params.sys.n*(params.ctrl.T+1)),solver_SDP);
% ctrl_aPWB.gamma-ctrl_aFBS.gamma
%% simulate the closed-loop system for nrSteps time steps and nrTraj
% Define noise simulations
% w_args{1} = {@zero, {params.sys.n, params.ctrl.T}};
% w_args{2} = {@constant, {1/sqrt(2),params.sys.n, params.ctrl.T}};
% w_args{3} = {@gaussian, {params.ctrl.T, zeros(params.sys.n,1), eye(params.sys.n)}};
% w_args{4} = {@uniform, {params.ctrl.T, zeros(params.sys.n,1),ones(params.sys.n,1)}};
% w_args{5} = {@sine, {params.sys.n, params.ctrl.T, 2/3*ones(params.sys.n), 0.1*ones(params.sys.n), zeros(params.sys.n)}};
% w_args{6} = {@sawt, {params.sys.n, params.ctrl.T, 0.5*ones(params.sys.n), 0.1*ones(params.sys.n), zeros(params.sys.n)}};
% 
% 
% %For all defined controllers and disturbance sims
% for i=1:length(ctrl)
%     for j=1:length(w_args)
%         [x(:,:,i,j), u(:,:,i,j), cost(:,i,j), cost_tot(:,i,j),w(:,:,i,j)] = sim_LQ(ctrl{i},w_args{j}{:},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
%     end
% end
% 
% %Compare all costs without NC
% for i=1:length(ctrl)
%     for j=1:length(w_args)
%         cost_comp(j,i)=cost_tot(end,i,j);
%     end
% end
% %Normalize
% for i=1:length(w_args)
%     cost_comp(i,:)=cost_comp(i,:)/min(cost_comp(i,:));
% end
% for i=1:length(ctrl)
%     ctrl_name{i}=ctrl{i}.name;
% end
% ctrl_name
% cost_comp

%%
if ctrl{1}.gamma-ctrl{3}.gamma>1e-5
    error('DR conservativeness reduced!')
end
if ctrl{2}.gamma-ctrl{4}.gamma>1e-5
    error('DR conservativeness reduced!')
end

%% plot results

% figIdx = 1:3;
% if params.plot.show
%     plot_x('state-time', figIdx(1), x_dr, params.sys.X, params.plot);
%     plot_x('state-state', figIdx(2), x_dr, params.sys.X, params.plot);
%     plot_u(figIdx(3), u_dr, sys.U, params.plot);
% end