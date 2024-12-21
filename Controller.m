classdef Controller < handle
    %CONTROLLER Abstract parent class for controller implementation
    %   Class to determine the input of a system given its current state
    
    properties
        params  %struct of parameters
        
    end
    
    methods
        function obj = Controller(sys, params)
            %CONTROLLER Construct an instance of this class
            %   Constructs a controller for a given system
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            % Initialize properties
            obj.params = params;
        end
    end
end