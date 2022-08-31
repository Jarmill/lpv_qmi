classdef lpvh2 < lpvstab
    %LPVSTAB stabilizing control of LPV systems using QMIs    
    
    %in experiments, stick with individual sample noise bounds.
    %until otherwise implemented
    
    properties
        
    end
    
    methods
        function obj = lpvh2(traj)
            %LPVSTAB Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@lpvstab(traj);         
        end
        
        function [out] = h2(obj, Th_vert, C, D, F)
            %h2: main h2 control routine and execution
            [cons, vars] = obj.make_program(Th_vert, C, D, F);
            out = obj.solve_program(cons, trace(vars.Z), vars);
        end
        
        %% form and solve the program
        function [cons, vars] = make_program(obj, Th_vert, C, D, F)
            %MAKE_PROGRAM form the LMI program in YALMIP
            
            %Th_vert: vertices of parameter region                        
                             
            n = obj.traj.n;
            Nv = size(Th_vert, 2);
            
            Nz = size(C, 1);
            
            Y = sdpvar(n, n);
            Z = sdpvar(Nz, Nz);

            %storage of variables
            vars = struct;
            vars.Y = Y;
            vars.M = cell(Nv, 1);
            vars.Cv = cell(Nv, 1);
            vars.Zmat = cell(Nv, 1);
            vars.a = zeros(Nv, 1, 'like', sdpvar);
            vars.b = zeros(Nv, 1, 'like', sdpvar);
            
            %specifically for h2 norm control
%             vars.gamma2 = sdpvar(1,1);
            vars.Z = Z;
            
            %storage of constraints
            cons = [(Y >= eye(obj.traj.n)*obj.delta); Z >= 0];
            %iterate through each subsystem and generate QMI
            for v = 1:Nv
                thv = Th_vert(:, v);
                
                [cons_vert, vars_vert] = obj.vertex_qmi(Y, Z, thv, C, D, F);
                
                cons = [cons; cons_vert];
                vars.a(v) = vars_vert.a;
                vars.b(v) = vars_vert.b;
                vars.M{v} = vars_vert.M;                
                vars.Cv{v} = vars_vert.Cv;
                vars.Zmat{v} = vars_vert.Zmat;
                
            end
            
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
            
            v2Y = kron(thv*thv', Y);
            vM = kron(thv, M');
            
            %form the controller matrix
            CT_bot = [-v2Y, -vM, zeros(L*n, n);
                      -vM', zeros(m, m), M;
                      zeros(n, L*n), M', Y];
                  
            CT = blkdiag(Y-b*eye(n)-F*F', CT_bot);
            Cv = CT - a*obj.Psi_E;
            
%             sCT = size(sCT, 1); %should be (L+2)n + m

            Zmat = [Z, C*Y + D*M; (C*Y + D*M)', Y];
            

            cons_vert = [Cv >= 0; a >= 0; b>=obj.delta; Zmat >= 0];
            vars_vert = struct('M', M, 'a', a, 'b', b, 'Cv', Cv, 'Zmat', Zmat);
        end
        
        function out = recover(obj, vars, sol)
            %RECOVER get the controllers and parameters
            out = struct;
            %variables
            out.sol = sol;
            out.a = value(vars.a);
            out.b = value(vars.b);
            
            %parameter-independent quadratic lyapunov function
            out.Y = value(vars.Y);
            
            Nv = length(vars.M);
            %controllers and qmi
            out.K = cell(Nv, 1);
            out.M = cell(Nv, 1);
            out.Cv = cell(Nv, 1);
            
            out.Z = value(vars.Z);
            out.h2 = sqrt(trace(out.Z));
%             out.gamma = sqrt(value(vars.gamma2));
            
            Nz = size(vars.Z, 1);
            out.eig_Y = eig(out.Y);
            sv = size(vars.Cv{1}, 1);
            out.eig_Cv = zeros(sv, Nv);
            out.eig_Zmat = zeros(size(vars.Zmat{1}, 1), Nv);
            for v = 1:Nv
                out.M{v} = value(vars.M{v});
                out.K{v} = (out.Y \ out.M{v}')';
                out.Cv{v} = value(vars.Cv{v});                
                out.eig_Cv(:, v) =eig(out.Cv{v});
                out.Zmat{v} = value(vars.Zmat{v});
                out.eig_Zmat(:, v) =eig(out.Zmat{v});
            end
        end
        
        function [valid, eig_out] = validate_stab_multi(obj, sys, Th_vert, K)
            if ~iscell(sys)
                sys = {sys};
            end
            Nsys = length(sys);
            valid = zeros(Nsys, 1);
            eig_out = cell(Nsys, 1);
            
            for i = 1:length(sys) 
                [valid(i), eig_out{i}] = obj.validate_stab(sys{i}, Th_vert, K);
            end
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
            
            valid = all(abs(eig_out) < 1, 'all');
        end
    end
end

