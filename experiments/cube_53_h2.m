
rng(45, 'twister');
n = 5;
m = 3;
L = 3;
% Tmax = 38; 
% Tmax = 36;
% Tmax = 50;
Tmax = 50;


LS = lpvsim(n, m, L);
cube = [1,1,1;
        1,1,-1;
        1,-1,1;
        1,-1,-1;
        -1,1,1;
        -1,1,-1;
        -1,-1,1;
        -1,-1,-1]';
        
Th_vert = diag([0.3, 0.3, 0.5])*cube + [0,0.5, 1]';     


% Th_vert = [1 0 0;
%            1 0.2 0;
%            1 0.3 0.5;
%            1.1 0 0;
%            1.1 0 0.5;
%            1.1 0.1 0.5]';

LS.sampler.th = @() sample_th(Th_vert);
traj = LS.sim(Tmax);


%run QMI solver
% LP = lpvstab(traj);
% 
% out = LP.stab(Th_vert);
% 
% out.sol.problem

Q = eye(n);
R = eye(m)*2;

C = [sqrt(Q); zeros(m, n)];
D = [zeros(n, m); sqrt(R)];
% F = [1; 0]*0.1;
F = eye(n)*0.2;

H2 = lpvh2(traj);


out = H2.h2(Th_vert, C, D, F);
disp(out.sol.problem)