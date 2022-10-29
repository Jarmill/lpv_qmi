%generate the trajectory
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

rng(40, 'twister');
n = 2;
m = 2;
L = 2;
Tmax = 35;

epsilon = 0.1;
% epsilon = 0.15;
LS = lpvsim(n, m, L, epsilon);
% traj = LS.sim(Tmax);

Th_vert = [-1 -1 1 1;
           1  -1 1 -1];

Th_shift = [1; 0];
Th_vert = Th_vert + Th_shift;
LS.sampler.th = @() (2*rand(2, 1) - 1) + Th_shift;
       
traj = LS.sim(Tmax);

%parameters 

Q = eye(n);
R = eye(m)*2;

C = [sqrt(Q); zeros(m, n)];
D = [zeros(n, m); sqrt(R)];
% F = [1; 0]*0.1;
F = eye(n);
       
% Th_vert = Th_vert + [0; 0.5];
%run QMI solver
if SOLVE
% LP = lpvstab(traj);
H2 = lpvh2(traj);

%constant controller
% H2.const_K = 1;

out = H2.h2(Th_vert, C, D, F);
disp(out.sol.problem)
end
