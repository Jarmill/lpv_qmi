%generate the trajectory
rng(45, 'twister');
n = 5;
m = 3;
L = 3;
% Tmax = 38; 
Tmax = 36;


LS = lpvsim(n, m, L);
% 
% Th_vert = [1 0 0;
%            1 0.2 0;
%            1 0.3 0.2;
%            1.1 0 0;
%            1.1 0 0.2;
%            1.1 0.1 0.2]';


Th_vert = [1 0 0;
           1 0.2 0;
           1 0.3 0.5;
           1.1 0 0;
           1.1 0 0.5;
           1.1 0.1 0.5]';

LS.sampler.th = @() sample_th(Th_vert);
traj = LS.sim(Tmax);


%run QMI solver
LP = lpvstab(traj);

out = LP.stab(Th_vert);

out.sol.problem

%% acquire samples

%% validate samples
if ~out.sol.problem
Nsys = 10;
sys_smp = LS.sample_sys(traj, Nsys);

    [valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end