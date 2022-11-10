classdef lpvh2_cont < lpvh2
    %LPVH2_CONT worst-case H2-optimal control of LPV systems using QMIs    
    %(continuous time)
    
    %in experiments, stick with individual sample noise bounds.
    %until otherwise implemented

    methods
        function obj = lpvh2_cont(traj)
            %LPVSTAB Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@lpvh2(traj);         
        end
                
        %% helper routines        
        function [cons_vert, vars_vert] = vertex_qmi(obj, Y, Z, thv,  C, D, F);
            %VERTEX_QMI: find the QMI associated with subsystem stabilization
            %at the parameter vertex thv            
            
            %find sizes
            n = obj.traj.n;
            m = obj.traj.m;
            L = size(thv, 1);
            
            %declare variables
            M = sdpvar(m, n, 'full'); %K = (Y \ M')'
            a = sdpvar(1,1);
            b = sdpvar(1,1);
            
            vY = kron(thv, Y);
            
             CT_row = -[vY; M];
%             v2Y = kron(thv*thv', Y);
%             vM = kron(thv, M');
            
            %form the controller matrix
            CT = [zeros(n)-b*eye(n)-F*F', CT_row'; CT_row, zeros(n*(L) + m)];
            
            Cv = CT - a*obj.Psi;
            
%             sCT = size(sCT, 1); %should be (L+2)n + m

            Zmat = [Z, C*Y + D*M; (C*Y + D*M)', Y];
            

            cons_vert = [Cv >= 0; a >= 0; b>=obj.delta; Zmat >= obj.delta];
            vars_vert = struct('M', M, 'a', a, 'b', b, 'Cv', Cv, 'Zmat', Zmat);
        end
        
    
       
        
      
    end
end

