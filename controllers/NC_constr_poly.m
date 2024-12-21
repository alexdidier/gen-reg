classdef NC_constr_poly
    %NC_CONSTR_POLY Optimal non-causal controller with constraints with
    %polytopic disturbance bounds
    %   Class for computing the optimal non-causal controller with state
    %   and input constraints for polytopic bounded disturbances
    
    properties
        name
        params
        sys
        M
        nc
        O
    end
    
    methods
        function obj = NC_constr_poly(sys, params, solver)
            %NC_CONSTR_POLY Construct an instance of this class
            %   Construct NC_constr_poly Class
            %%% Parse inputs %%%
            switch nargin
                case 2
                    solver = '';
                    
                case 3
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            obj.name='NC_constr_poly';
            obj.nc = true;
            obj.params = params;
            obj.sys = sys;
            
            % Some definitions for readability
            n = obj.sys.n;
            m = obj.sys.m;
            A = obj.sys.A;
            B = obj.sys.B;
            Q = obj.params.Q;
            R = obj.params.R;
            T = obj.params.T;
            x_0 = obj.params.x_0;
            
            % System matrices used in the SLS formulations
            ZA = kron(diag(ones(1,T),-1),A);
            ZB = kron(diag(ones(1,T),-1),B);
            Q_cal = kron(eye(T+1),Q);
            R_cal = kron(eye(T+1),R);
            C=[Q_cal, zeros(size(Q_cal,1),size(R_cal,2));zeros(size(Q_cal,1),size(R_cal,2))', R_cal ];
            
            %State and input constraints
            H_x = obj.sys.H_x;
            h_x = obj.sys.h_x;
            H_u = obj.sys.H_u;
            h_u = obj.sys.h_u;
            
            H_z = [kron(eye(T+1),H_x),zeros(size(H_x,1)*(T+1),m*(T+1)); zeros(size(H_u,1)*(T+1),n*(T+1)),kron(eye(T+1),H_u)];
            h_z = [kron(ones(T+1,1),h_x);kron(ones(T+1,1),h_u)];
            
            H_w = [kron(eye(T),obj.sys.H_w)];
            h_w = [kron(ones(T,1),obj.sys.h_w)];
            % Optimisation variables
            
            Phi_x = sdpvar(n*(T+1),n*(T+1),'full');
            Phi_u = sdpvar(m*(T+1), n*(T+1),'full');
            Gamma = sdpvar(size(H_z,1),size(H_w,1),'full');
            
            % Constraints
            
            constraints = [Gamma(:)>=0];

            Phi = [Phi_x; Phi_u];
            Phi_0 = Phi(:,1:n);
            Phi_w = Phi(:,n+1:end);
            
            constraints = [constraints, [eye(size(ZA))-ZA, -ZB]*[Phi_x; Phi_u] == eye(size(ZA)) ];

            constraints = [constraints, Gamma*h_w<=h_z-H_z*Phi_0*x_0, Gamma*H_w==H_z*Phi_w];
            
            % Objective
            
            objective = norm(C^0.5*Phi,'fro');
            
            % Solve the SDP
            
            info=optimize(constraints,objective,sdpsettings('solver',solver,'verbose',0,'savesolveroutput',1));
            
            if (info.problem ~=0)
                error(info.info)
            end
            
            obj.M = value(Phi_u);    
            obj.O = value(Phi)'*C*value(Phi);
        end
        
        function u = inp(obj, w, k)
            %INP Returns input to be applied
            % Computes input based on controller matrix K and current time
            % step
            
            %%% Parse inputs %%%
            switch nargin
                case 3
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            m = obj.sys.m;
            u = obj.M((k-1)*m+1:k*m,:)*vec(w);
            
        end
    end
end