classdef lpvsim
    %LPVSIM_OPTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n = 3;      %number of states
        m = 2;      %number of inputs
        L = 2;      %number of parameters

        sampler = [];

        A_scale = 1.2;      %perturb randomly generated stable ss systems to make them possibly unstable
    end

    methods

        function obj = lpvsim(n, m, L)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            if nargin == 3
                obj.n = n;
                obj.m = m;
                obj.L = L;
            end

            obj.sampler = struct('th', @() 2*rand(obj.L,1)-1, ...
                                 'u', @() 2*rand(obj.m, 1)-1, ...
                                 'w', @() normalize(randn(obj.n,1), 1, 'norm', 2));

            %             obj.Property1 = inputArg1 + inputArg2;
        end

        function out = sim(obj, T, sys, epsilon, x0)
            %simulate an LPV trajectory with individual sample noise bound
            %epsilon (L2 norm)
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

            if nargin < 5
                epsilon = 0.1;
            end

            
            X = [x0, zeros(obj.n, T)];
            U = zeros(obj.m, T);
            W_true = zeros(obj.n, T);
            Th = zeros(obj.L, T);
            %main simulation loop
            xcurr = x0;
            for i = 1:T
                %inputs
                ucurr = obj.sampler.u();
                wcurr = obj.sampler.w()*epsilon;
                thcurr = obj.sampler.th();
                
                %propagation
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
            out.epsilon = epsilon;            
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.L = obj.L;

        end

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


    
    end
%     function sa

end

