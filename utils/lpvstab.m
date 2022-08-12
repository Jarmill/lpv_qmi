classdef lpvstab
    %LPVSTAB stabilizing control of LPV systems using QMIs    
    
    %in experiments, stick with individual sample noise bounds.
    %until otherwise implemented
    
    properties
        
        %start out with only one trajectory, then expand to multiple
        out; %trajectory trace
        
        Psi; %data matrix
        
        vars;
        
        delta = 1e-5; %tolerance for being positive (positive definite)
                
    end
    
    methods
        function obj = lpvstab(out)
            %LPVSTAB Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj.out = out;
            
            obj.Psi = sample_matrix(out.X, out.U, out.epsilon, out.Th);
            
%             vars = [];
        end
        
        function [cons, vars] = make_program(obj, Th_vert)
            %MAKE_PROGRAM form the LMI program in YALMIP
            
            %Th_vert: vertices of parameter region
            
            [L, Nv] = size(Th_vert);

            [cons_vars, vars] = obj.make_vars(Nv);
            
        end
        
%         function [cons_vars, vars] = obj.make_vars(obj, Nv)
%             %declare all variab
%             n = size(obj.out.X, 1);
%             m = size(obj.out.U, 1);
%             L = size(obj.out.L, 1);
%             
%             
%             
%         end
        
    end
end

