classdef GenReg_Poly
    %GENREG_PWB Generalised dynamic regret controller for polytopic PWB dist.
    %   Class for solving the generalised dynamic regret problem for
    %   pointwise in time bounded polytopic disturbances
    
    properties
        name
        params
        sys
        K
        gamma
        lambda
    end
    
    methods
        function obj = GenReg_Poly(sys, params, weight, solver)
            %GENREG_PWB Construct an instance of this class
            %   Construct GenReg Class
            %%% Parse inputs %%%
            switch nargin
                case 2
                    weight='dr';
                    solver = '';
                    
                case 3
                    solver = '';
                    
                case 4
                    
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            
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
            omega = obj.sys.omega;
            
            % System matrices used in the SLS formulations
            
            ZA = kron(diag(ones(1,T),-1),A);
            ZB = kron(diag(ones(1,T),-1),B);
            Q_cal = kron(eye(T+1),Q);
            R_cal = kron(eye(T+1),R);
            C=[Q_cal, zeros(size(Q_cal,1),size(R_cal,2));zeros(size(Q_cal,1),size(R_cal,2))', R_cal ];
            C_sqrt=C^0.5;
            F=zeros((T+1)*n,(T+1)*m);
            G=zeros((T+1)*n,(T+1)*n);
            
            for i=1:T
                for j=1:T+1
                    if j<=i
                        F(n*i+1:n*(i+1),m*(j-1)+1:m*j)=A^(i-j)*B;
                        G(n*i+1:n*(i+1),n*(j)+1:n*(j+1))=A^(i-j)*eye(n);
                    end
                end
                G((i-1)*n+1:i*n,1:n)=A^(i-1);
            end
            G(T*n+1:(T+1)*n, 1:n) = A^T;
            O=G'/(inv(Q_cal)+F/R_cal*F')*G;
            O_1 = O(1:n,1:n);
            O_2 = O(n+1:end, 1:n);
            O_3 = O(n+1:end, n+1:end);
            
            if ischar(weight)
                switch weight
                    case 'dr'
                        obj.name='DR_PWB';
                        W_w=eye(n*(T+1));
                    case 'cr'
                        obj.name='CR_PWB';
                        W_w=O;
                    otherwise
                        error('Not a valid generalised regret case.')
                end
            else
                obj.name='GR_PWB';
                W_w=weight;
            end
            assert(all(eig(W_w)>0),'Disturbance weight is not positive definite.');
            
            W_1 = W_w(1:n,1:n);
            W_2 = W_w(n+1:end, 1:n);
            W_3 = W_w(n+1:end, n+1:end);
            
            % Optimisation variables
            
            gamma = sdpvar(1,1);
            lambda = sdpvar(T*size(obj.sys.H_w,1),1);
            lambda_w=sdpvar(1,1);
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
            
            lambdah=zeros(1,1);
            for i=1:T
                for j=1:size(obj.sys.h_w,1)
                    lambdah = lambdah + lambda((i-1)*size(obj.sys.h_w,1)+j)*obj.sys.h_w(j);
                end
            end
            lambdaH = zeros(1,n*T);
            for i=1:T
                for j=1:size(obj.sys.H_w,1)
                temp=zeros(1,n*T);
                temp(1,(i-1)*n+1:i*n)=obj.sys.H_w(j,:);
                lambdaH = lambdaH + lambda((i-1)*size(obj.sys.h_w,1)+j)*temp;
                end
            end
            
            constraints = [constraints, lambda(:)>=0,lambda_w>=0];
            constraints = [constraints, [eye(size(ZA))-ZA, -ZB]*[Phi_x; Phi_u] == eye(size(ZA)) ];
            constraints = [constraints, [ x_0'*(O_1+gamma*W_1)*x_0-lambdah-lambda_w*omega, x_0'*(O_2'+gamma*W_2')+0.5*lambdaH, x_0'*Phi_0'*C_sqrt; ...
                (O_2+gamma*W_2)*x_0+0.5*lambdaH', lambda_w*eye(size(O_3))+O_3+gamma*W_3, Phi_w'*C_sqrt; ...
                C_sqrt'*Phi_0*x_0, C_sqrt'*Phi_w, eye(size(C))]>=0];

            % Objective
            
            objective = gamma;
            
            % Solve the SDP
            
            info=optimize(constraints,objective,sdpsettings('solver',solver,'verbose',0,'savesolveroutput',1));
            
            if (info.problem ~=0)
                error(info.info)
            end
            
            obj.gamma=value(gamma);
            obj.lambda=value(lambda);
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