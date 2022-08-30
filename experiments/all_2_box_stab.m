%generate the trajectory
SOLVE = 1;
SAMPLE = 0;
PLOT = 1;

rng(40, 'twister');
n = 2;
m = 2;
L = 2;
Tmax = 35;
% Tmax = 45;
% 

% n = 3;
% m = 2;
% L = 2;
% Tmax = 35;


epsilon = 0.1;
% epsilon = 0.25;
LS = lpvsim(n, m, L, epsilon);
% traj = LS.sim(Tmax);
traj = LS.sim(Tmax);

Th_vert = [-1 -1 1 1;
           1  -1 1 -1];

Th_vert = Th_vert + [0; 0.5];
%run QMI solver
if SOLVE
LP = lpvstab(traj);

out = LP.stab(Th_vert);
disp(out.sol.problem)
end
if SAMPLE
%% acquire samples
Nsys = 10;
sys_smp = LS.sample_sys(traj, Nsys);

%% validate samples
if ~isempty(out)
[valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end
end

%% plot sample trajectories
if PLOT
    %% test out the interpolating controller
[L, Nv] = size(Th_vert);


% th = sdpvar(L, 1);
% c = sdpvar(Nv, 1);
% 
% cons = [c>=0; sum(c)==1; Th_vert*c==th];
% Pth = optimizer(cons, 0, sdpsettings('solver', 'mosek'), th, c);
Pth = get_vertex_interp(Th_vert);
th_test = [1; 0.4];

Kth = @(th) K_interp(Pth, out.K, th);

% uth = @(th, x) Kth(th)*x;


LS.sampler.u = @(th,x) Kth(th)*x;

%% plot a trajectory

Ntraj = 100;

x0 = [2; 1];
figure(3)
clf
hold on
for i = 1:Ntraj
    traj_closed = LS.sim(Tmax, traj.ground_truth, x0);

    plot(traj_closed.X(1, :), traj_closed.X(2, :))
end
scatter(traj_closed.X(1, 1), traj_closed.X(2,1), 200, 'ok')

xlabel('x_1')
ylabel('y_1')


end