%generate the trajectory
rng(3, 'twister');
n = 2;
m = 2;
L = 2;
LS = lpvsim(n, m, L);
traj = LS.sim(35);

%define vertices
% Th_vert = [-1,1];
Th_vert = [-1 -1 1 1;
           1  -1 1 -1];


%run QMI solver
LP = lpvstab(traj);

out = LP.stab(Th_vert);