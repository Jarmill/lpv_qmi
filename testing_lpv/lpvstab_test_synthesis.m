%generate the trajectory
rng(40, 'twister');
n = 2;
m = 2;
L = 2;
LS = lpvsim(n, m, L);
Tmax = 35;
traj = LS.sim(Tmax);

%define vertices
% Th_vert = [-1,1];
Th_vert = [-1 -1 1 1;
           1  -1 1 -1];


%run QMI solver
LP = lpvstab(traj);

out = LP.stab(Th_vert);

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
traj_closed = LS.sim(Tmax, traj.ground_truth, 0, x0);

plot(traj_closed.X(1, :), traj_closed.X(2, :))
end
scatter(traj_closed.X(1, 1), traj_closed.X(2,1), 200, 'ok')

xlabel('x_1')
ylabel('y_1')
% K_test = K_interp(Pth, out.K, th_test);

% c_test = Pth(th_test);
% 
% K_test = zeros(size(out.K{1}));
% for i = 1:Nv
%     K_test = out.K{i}*c_test(i);
% end
% cons = [c>=0; sum(c)==1; Th_vert*c==th_test];
% sol_th = optimize(cons, 0, []);
% th_rec_test = 