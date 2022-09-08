classdef lpvstab_cont < lpvstab
    %LPVSTAB stabilizing control of LPV systems using QMIs    
    
    %in experiments, stick with individual sample noise bounds.
    %until otherwise implemented
    
    properties
                
    end
    
    methods
        function obj = lpvstab_cont(traj)
            %LPVSTAB Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            obj@lpvstab(traj);
            
        end
        
        %% helper routines        
        function [cons_vert, vars_vert] = vertex_qmi(obj, Y, thv);
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
%             CT_bot = [-v2Y, -vM, zeros(L*n, n);
%                       -vM', zeros(m, m), M;
%                       zeros(n, L*n), M', Y];
%                   
%             CT = blkdiag(Y-b*eye(n), CT_bot);
            CT = [zeros(n), CT_row'; CT_row, zeros(n*(L) + m)];

            Cv = CT - a*obj.Psi;
            
%             sCT = size(sCT, 1); %should be (L+2)n + m
            cons_vert = [Cv >= 0; a >= 0; b>=obj.delta];
            vars_vert = struct('M', M, 'a', a, 'b', b, 'Cv', Cv);
        end
               
        function [valid, eig_out] = validate_stab(obj, sys, Th_vert, K)
            %validate stability of every vertex given stab and the
            %controller K
            
            n = size(sys.A{1}, 1);
            [L, Nv] = size(Th_vert);
            eig_out = zeros(n, Nv);
            
            for v = 1:Nv
                Acl = sys.B * K{v};
                for i = 1:L
                    Acl = Acl + sys.A{i}*Th_vert(i, v);
                end
                eig_out(:, v) = eig(Acl);
            end
            
            valid = all(real(eig_out) < 0, 'all');
        end
    end
end

