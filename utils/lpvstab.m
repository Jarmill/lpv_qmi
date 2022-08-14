classdef lpvstab
    %LPVSTAB stabilizing control of LPV systems using QMIs    
    
    %in experiments, stick with individual sample noise bounds.
    %until otherwise implemented
    
    properties
        
        %start traj with only one trajectory, then expand to multiple
        traj; %trajectory trace
        
        Psi; %data matrix 
        Psi_E; %data matrix with zero padding
        
%         vars;
        
        delta = 1e-5; %tolerance for being positive (positive definite)
       
        
        opts = sdpsettings('solver', 'mosek');
    end
    
    methods
        function obj = lpvstab(traj)
            %LPVSTAB Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj.traj = traj;
            
            obj.Psi = sample_matrix(traj.X, traj.U, traj.epsilon, traj.Th);
            
            obj.Psi_E = blkdiag(obj.Psi, zeros(obj.traj.n));            
        end
        
        function [out] = stab(obj, Th_vert)
            %STAB: main stabilization routine and execution
            [cons, vars] = obj.make_program(Th_vert);
            out = obj.solve_program(cons, vars);
        end
        
        %% form and solve the program
        function [cons, vars] = make_program(obj, Th_vert)
            %MAKE_PROGRAM form the LMI program in YALMIP
            
            %Th_vert: vertices of parameter region
                             
            n = obj.traj.n;
            Nv = size(Th_vert, 2);
            
            Y = sdpvar(n, n);

            %storage of variables
            vars = struct;
            vars.Y = Y;
            vars.M = cell(Nv, 1);
            vars.Cv = cell(Nv, 1);
            vars.a = zeros(Nv, 1, 'like', sdpvar);
            vars.b = zeros(Nv, 1, 'like', sdpvar);
            
            %storage of constraints
            cons = (Y >= eye(obj.traj.n)*obj.delta);
            %iterate through each subsystem and generate QMI
            for v = 1:Nv
                thv = Th_vert(:, v);
                
                [cons_vert, vars_vert] = obj.vertex_qmi(Y, thv);
                
                cons = [cons; cons_vert];
                vars.a(v) = vars_vert.a;
                vars.b(v) = vars_vert.b;
                vars.M{v} = vars_vert.M;                
                vars.Cv{v} = vars_vert.Cv;
                
            end
            
        end
        
        function [out] = solve_program(obj, cons, vars)
            %run the program
            sol = optimize(cons, 0, obj.opts);
            out.sol = sol;
            
            if sol.problem==0
                %successful trajectory execution
                out = obj.recover(vars, sol);
            end
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
            
            v2Y = kron(thv*thv', Y);
            vM = kron(thv, M');
            
            %form the controller matrix
            CT_bot = [-v2Y, -vM, zeros(L*n, n);
                      -vM', zeros(m, m), M;
                      zeros(n, L*n), M', Y];
                  
            CT = blkdiag(Y-b*eye(n), CT_bot);
            Cv = CT - a*obj.Psi_E;
            
%             sCT = size(sCT, 1); %should be (L+2)n + m
            cons_vert = [Cv >= 0; a >= 0; b>=obj.delta];
            vars_vert = struct('M', M, 'a', a, 'b', b, 'Cv', Cv);
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
            
            
            out.eig_Y = eig(out.Y);
            sv = size(vars.Cv{1}, 1);
            out.eig_Cv = zeros(sv, Nv);
            for v = 1:Nv
                out.M{v} = value(vars.M{v});
                out.K{v} = (out.Y \ out.M{v}')';
                out.Cv{v} = value(vars.Cv{v});                
                out.eig_Cv(:, v) =eig(out.Cv{v});
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

