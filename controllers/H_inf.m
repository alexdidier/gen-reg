classdef H_inf
    %H_inf Optimal H_inf controller
    %   Class for solving the H_inf problem
    
    properties
        name
        params
        sys
        K
    end
    
    methods
        function obj = H_inf(sys, params, solver)
            %H_inf Construct an instance of this class
            %   Construct H_inf Controller Class
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    solver = '';
                    
                case 3
                    
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            obj.name='Hinf';
            obj.params = params;
            obj.sys = sys;
            
            % Some definitions for readability
            n = obj.sys.n;
            m = obj.sys.m;
            A = obj.sys.A;
            B = obj.sys.B;
            omega = obj.sys.omega;
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
            Cinv=inv(C);
            
            % Optimisation variables
            
            gamma = sdpvar(1,1);
            lambda = sdpvar(1,1);
            Phi_x = sdpvar(n*(T+1),n*(T+1),'full');
            Phi_u = sdpvar(m*(T+1), n*(T+1),'full');
            
            % Constraints
            
            constraints = [];
            
            for i=1:T
                for j=i+1:T+1
                    Phi_x((i-1)*n+1:i*n,(j-1)*n+1:j*n)=zeros(n,n);
                    Phi_u((i-1)*m+1:i*m,(j-1)*n+1:j*n)=zeros(m,n);
                end
            end
            Phi_u(T*m+1:(T+1)*m,:)=zeros(size(Phi_u(T*m+1:(T+1)*m,:)));
            Phi = [Phi_x;Phi_u];
            Phi_0 = Phi(:,1:n);
            Phi_w = Phi(:,n+1:end);
            
            
            constraints = [constraints, lambda>=0];
            constraints = [constraints, [eye(size(ZA))-ZA, -ZB]*Phi == eye(size(ZA)) ];
            constraints = [constraints, [ gamma-lambda*omega, zeros(1,size(Phi_w,2)), x_0'*Phi_0';
                              zeros(size(Phi_w,2),1), lambda*eye(size(Phi_w,2)), Phi_w'; ...
                              Phi_0*x_0, Phi_w, Cinv]>=0];

            % Objective
            
            objective = gamma;
            
            % Solve the SDP
            
            info=optimize(constraints,objective,sdpsettings('solver',solver,'verbose',0,'savesolveroutput',1));
            
            if (info.problem ~=0)
                error(info.info)
            end
            
            obj.K = value(Phi_u)/value(Phi_x);
        end
        
        function u = inp(obj, x, k)
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
            u = obj.K((k-1)*m+1:k*m,:)*vec(x);
            
        end
    end
end