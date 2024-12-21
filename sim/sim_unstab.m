%% Setup

clc;
clear all;
close all;

% fix random number generator
rng(1);

% Select SDP Solver
solver_SDP = 'mosek';

%% get parameters and define system and controller

params = get_params('params_unstab');

ctrl{1} = NC(params.sys, params.ctrl);
ctrl{2} = H_2(params.sys, params.ctrl, solver_SDP);
ctrl{3} = H_inf(params.sys,params.ctrl, solver_SDP);
tic
ctrl{4} = GenReg(params.sys, params.ctrl,'dr', solver_SDP);
toc
%ctrl{5} = GenReg(params.sys, params.ctrl, 'cr', solver_SDP);

%ctrl{6} = GenReg_PWB(params.sys,params.ctrl,'dr',solver_SDP);
%ctrl{7} = GenReg_PWB(params.sys,params.ctrl,'cr',solver_SDP);

%% Controllers with state and input constraints

% ctrl_s{1} = NC_constr_poly(params.sys, params.ctrl, solver_SDP);
% ctrl_s{2} = H_2_constr_poly(params.sys, params.ctrl, solver_SDP);
% ctrl_s{3} = H_inf_constr_poly(params.sys,params.ctrl, solver_SDP);
% 
% params.ctrl.O = ctrl_s{1}.O;
% ctrl_s{4} = GenReg_constr_poly(params.sys,params.ctrl,'dr',solver_SDP);
% ctrl_s{5} = GenReg_constr_poly(params.sys,params.ctrl,'cr',solver_SDP);
% %params.ctrl.Phi_BE=ctrl_s{4}.Phi;
% ctrl_s{6} = GenReg_PWB_constr_poly(params.sys,params.ctrl,'dr',solver_SDP);
% ctrl_s{7} = GenReg_PWB_constr_poly(params.sys,params.ctrl,'cr',solver_SDP);

%% simulate the closed-loop system for nrSteps time steps and nrTraj
% Define noise simulations
%w_args{1} = {@gaussian_trunc, {params.ctrl.T, zeros(params.sys.n,1), eye(params.sys.n),-ones(params.sys.n,1),ones(params.sys.n,1)},true};
%w_args{2} = {@uniform, {params.ctrl.T, zeros(params.sys.n,1),ones(params.sys.n,1)},true};
%w_args{3} = {@uniform, {params.ctrl.T, 0.5*ones(params.sys.n,1),ones(params.sys.n,1)},true};
w_args{1} = {@constant, {1,params.sys.n, params.ctrl.T},false};
%w_args{5} = {@sine, {params.sys.n, params.ctrl.T, ones(params.sys.n), 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
%w_args{6} = {@sawt, {params.sys.n, params.ctrl.T, ones(params.sys.n), 0.1*ones(params.sys.n), zeros(params.sys.n)},false};
%w_args{7} = {@step, {params.sys.n, params.ctrl.T}, false};
%w_args{8} = {@stair, {params.sys.n, params.ctrl.T}, false};
%w_args{8} = {@adversarial, {ctrl_s{3}.Phi, ctrl_s{3}.C, params.sys.n, params.ctrl.T, params.sim.x_0, 'Poly', params.sys.H_w}}

%For all defined controllers and disturbance sims
for i=1:length(ctrl)
    for j=1:length(w_args)
        if w_args{j}{3}==true
            %Ensure same disturbances for every different controller
            rng(1)
            for k=1:params.sim.nrTraj
                [x(:,:,i,j,k), u(:,:,i,j,k), cost(:,i,j,k), cost_tot(:,i,j,k),w(:,:,i,j,k)] = sim_LQ(ctrl{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
            end
            cost_tot(end,i,j,1)=sum(cost_tot(end,i,j,:))/params.sim.nrTraj;
        else 
            [x(:,:,i,j,1), u(:,:,i,j,1), cost(:,i,j,1), cost_tot(:,i,j,1),w(:,:,i,j,1)] = sim_LQ(ctrl{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
        end
    end
end

% %For all defined safe controllers and disturbance sims
% for i=1:length(ctrl_s)
%     for j=1:length(w_args)
%          if w_args{j}{3}==true
%             %Ensure same disturbances for every different controller
%             rng(1)
%             for k=1:params.sim.nrTraj
%                 [x_s(:,:,i,j,k), u_s(:,:,i,j,k), cost_s(:,i,j,k), cost_tot_s(:,i,j,k+1),w_s(:,:,i,j,k)] = sim_LQ(ctrl_s{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
%             end
%             cost_tot_s(end,i,j,1)=sum(cost_tot_s(end,i,j,2:end))/params.sim.nrTraj;
%         else 
%             [x_s(:,:,i,j,1), u_s(:,:,i,j,1), cost_s(:,i,j,1), cost_tot_s(:,i,j,1),w_s(:,:,i,j,1)] = sim_LQ(ctrl_s{i},w_args{j}{1:2},params.sys.A,params.sys.B,params.ctrl.Q,params.ctrl.R,params.sim.x_0,params.sim.nrSteps,false);
%         end
%     end
% end

%% Compare all costs without NC
for i=1:length(ctrl)-1
    for j=1:length(w_args)
        cost_comp(j,i)=cost_tot(end,i+1,j,1);
    end
end
%Normalize
% for i=1:length(w_args)
%     cost_comp(i,:)=cost_comp(i,:)/min(cost_comp(i,:));
% end
for i=1:length(ctrl)-1
    ctrl_name{i}=ctrl{i+1}.name;
end
ctrl_name
cost_comp

% %% Safe controllers
% for i=1:length(ctrl_s)-1
%     for j=1:length(w_args)
%         cost_comp_s(j,i)=cost_tot_s(end,i+1,j,1);
%     end
% end
% 
% %Normalize
% for i=1:length(w_args)
%     cost_comp_s(i,:)=cost_comp_s(i,:)/min(cost_comp_s(i,:));
% end
% 
% for i=1:length(ctrl_s)-1
%     ctrl_name_s{i}=ctrl_s{i+1}.name;
% end
% 
% ctrl_name_s
% cost_comp_s

%% plot results

% figIdx = 1:3;
% if params.plot.show
%     plot_x('state-time', figIdx(1), x_dr, params.sys.X, params.plot);
%     plot_x('state-state', figIdx(2), x_dr, params.sys.X, params.plot);
%     plot_u(figIdx(3), u_dr, sys.U, params.plot);
% end