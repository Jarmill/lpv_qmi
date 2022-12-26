classdef lpvsim
    %LPVSIM_OPTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n = 3;      %number of states
        m = 2;      %number of inputs
        L = 2;      %number of parameters
        epsilon = 0.1;

        sampler = [];

        A_scale = 1.2;      %perturb randomly generated stable ss systems to make them possibly unstable
    end

    methods

        function obj = lpvsim(n, m, L, epsilon)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            if nargin >= 3
                obj.n = n;
                obj.m = m;
                obj.L = L;
            end
            
            if nargin ==4
                obj.epsilon = epsilon;
            end
                

            obj.sampler = struct('th', @() 2*rand(obj.L,1)-1, ...
                                 'u', @(th, x) 2*rand(obj.m, 1)-1, ...
                                 'w', @() ball_sample(1, obj.n)');
%                                  'w', @() normalize(randn(obj.n,1), 1, 'norm', 2));

            %             obj.Property1 = inputArg1 + inputArg2;
        end

        function out = sim(obj, T, sys, x0)
            %simulate an continuous-time LPV trajectory with individual 
            %sample noise bound epsilon (L2 norm)
            out = struct;

            if nargin < 2
                T = 15;
            end

            if nargin < 3
                sys = obj.rand_sys();
            end
            if nargin < 4
                x0 = zeros(obj.n, 1);
                x0(1) = 1;
            end           

            
            X = [x0, zeros(obj.n, T)];
            U = zeros(obj.m, T);
            W_true = zeros(obj.n, T);
            Th = zeros(obj.L, T);
            %main simulation loop
            xcurr = x0;
            for i = 1:T
                %inputs
                
                wcurr = obj.sampler.w()*obj.epsilon;
                thcurr = obj.sampler.th();
                ucurr = obj.sampler.u(thcurr, xcurr);
                
                %propagation 
                %this is where the noise enters (process ?)
                xnext = sys.B*ucurr + wcurr;
                for k = 1:obj.L
                    xnext = xnext + sys.A{k}*xcurr*thcurr(k);
                end                
                
                %storage
                X(:, i+1) = xnext;
                U(:, i) = ucurr;
                W_true(:, i) = wcurr;
                Th(:, i) = thcurr;
                xcurr = xnext;
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            ground_truth.W = W_true;
            ground_truth.W2 = sum(W_true.^2, 1);
            
%             struct('W', W_true, 'A', sys.A, 'B', sys.B)

            %package up the output
            out.X = X;
            out.U = U;
            out.Th = Th;
            out.epsilon = obj.epsilon;            
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.L = obj.L;

        end
        
        function out = sim_closed_cont(obj, sys, Kth, x0, T, mu)
            %simulate an continuous-time controlled LPV trajectory 
            %there is no sample noise here.
            
            %Kth(th): gain-scheduled controller at parameter value theta
            out = struct;

            if nargin < 5
                T = 15;
            end
            if nargin < 6
                mu = 0.3;
            end

            if nargin < 2
                sys = obj.rand_sys();
            end
            
            if nargin < 3
                Kth = @(th) 0;
            end
            
            if nargin < 4
                x0 = zeros(obj.n, 1);
                x0(1) = 1;
            end           

            X = [];
            Th = [];
            U = [];
            Tlog = [];
            
%             X = [x0, zeros(obj.n, T)];
%             U = zeros(obj.m, T);
%             W_true = zeros(obj.n, T);
%             Th = zeros(obj.L, T);
%             %main simulation loop
%             traj = {};
            xprev = x0;
%             trajcurr = struct;
            tmax_curr = exprnd(mu);          
            t_all = 0;
            switch_times = 0;
            while t_all < T
                tmax_curr = exprnd(mu);
                tmax_curr = min(t_all + tmax_curr, T) - t_all;
                %inputs
                
%                 wcurr = obj.sampler.w()*obj.epsilon;
                thcurr = obj.sampler.th();
                Kcurr = Kth(thcurr);
                
                %propagation 
                Ath = 0;
                for ell = 1:length(thcurr)
                    Ath = sys.A{ell}*thcurr(ell);
                end
                
                Acl = Ath + sys.B*Kcurr;
                [tcurr, xcurr] = ode45(@(t, x) Acl*x, [0, tmax_curr], xprev);
                

     
                X = [X, xcurr'];
                U = [U, Kcurr*xcurr'];
                Th = [Th, thcurr];
                Tlog = [Tlog; t_all + tcurr];
                %storage
%                 X(:, i+1) = xnext;
%                 U(:, i) = ucurr;
%                 W_true(:, i) = wcurr;
%                 Th(:, i) = thcurr;
                xprev = xcurr(end, :)';
                switch_times = [switch_times; t_all + tmax_curr];
                t_all = t_all + tmax_curr;
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            
%             struct('W', W_true, 'A', sys.A, 'B', sys.B)

            %package up the output
            out.X = X;
            out.U = U;
            out.Th = Th;    
            out.t = Tlog;
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.L = obj.L;
            out.switch_times  = switch_times ;

        end

        %% generate sample plants inside the consistency set
        
        function sys_lpv = rand_sys(obj, A_scale)

            %randomly generate the LPV system
            if nargin < 2
                A_scale = obj.A_scale;
            end
            A = cell(obj.L, 1);
            for i = 1:obj.L
                sys= drss(obj.n, obj.n, obj.m);
                A{i} = A_scale*sys.a; %raise the likelihood that some systems are unstable
                if i == 1
                    B = sys.b;
                end
            end
            
            %package the output

            sys_lpv = struct;
            sys_lpv.A = A;
            sys_lpv.B = B;
        end
        
        function sys_smp = sample_sys(obj, traj, Nsys)
            %sample a set of systems consistent with the data
            %optima of linear objective over the quadratic region
            if nargin <3
                Nsys = 1;
            end
            n = size(traj.X,1);
            [m, T] = size(traj.U);
            L = size(traj.Th, 1);
            
            Xp = traj.X(:, 2:end);
            Xn = traj.X(:, 1:end-1);
            
            %declare variables and error object
            B = sdpvar(n, m, 'full');
            A = cell(L, 1);
            W = Xp -B*traj.U;
            for k = 1:obj.L
                A{k} = sdpvar(n, n, 'full');
                W = W - A{k}*(traj.Th(k, :) .* Xn);
            end                
            
            cons = [];
            for t = 1:T
                cons = [cons; norm(W(:, t), 2) <= traj.epsilon]; 
            end
%             W2 = sum(W.^2, 1);

%             cons = (W2' <= traj.epsilon^2);
            
            sys_smp = cell(Nsys, 1);
            %I'm not sure why the yalmip optimizer isn't working here.
%             opts = sdpsettings('solver','mosek', 'verbose', 2);
            opts = sdpsettings('solver','mosek', 'verbose', 0);

            for i = 1:Nsys
%                 sys_curr = struct('A', [], 'B', []);
                CB = randn(n, m);
                
                objective = sum(CB.*B, 'all');
                for k = 1:obj.L
                    CA = randn(n, n);
                    objective = objective + sum(CA.*A{k}, 'all');
                end
                %solve the program                
                sol = optimize(cons, objective, opts);
                
                %recover the solution
                sys_curr = struct;
                sys_curr.B = value(B);
                sys_curr.A = cellfun(@value, A, 'UniformOutput', false);
%                 sys_curr.W2 = value(W2);
                sys_curr.W = value(W);
                sys_curr.Wnorm = sqrt(sum(value(W).^2, 1));
                sys_smp{i} = sys_curr;
            end           
        end
              


    
    end
%     function sa

end

