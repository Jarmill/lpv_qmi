%generate the trajectory
rng(40, 'twister');
n = 2;
m = 2;
L = 2;
Tmax = 35;
% 

% n = 3;
% m = 2;
% L = 2;
% Tmax = 35;



LS = lpvsim(n, m, L);
traj = LS.sim(Tmax);

Th_vert = [-1 -1 1 1;
           1  -1 1 -1];


%run QMI solver
LP = lpvstab(traj);

out = LP.stab(Th_vert);


%% acquire samples
Nsys = 10;
sys_smp = LS.sample_sys(traj, Nsys);

%% validate samples
if ~isempty(out)
[valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end