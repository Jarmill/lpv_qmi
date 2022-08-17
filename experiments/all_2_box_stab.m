%generate the trajectory
SOLVE = 1;
SAMPLE = 0;

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