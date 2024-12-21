classdef NC
    %NC Optimal non-causal controller
    %   Class for computing the optimal non-causal controller
    
    properties
        name
        params
        sys
        M
        nc
    end
    
    methods
        function obj = NC(sys, params)
            %NC Construct an instance of this class
            %   Construct NC Class
            %%% Parse inputs %%%
            switch nargin
                case 2

                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            obj.name='NC';
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
            
            % System matrices
            Q_cal = kron(eye(T+1),Q);
            R_cal = kron(eye(T+1),R);

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
            obj.M = -(R_cal+F'*Q_cal*F)^(-1)*F'*Q_cal*G;     
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