%generate the trajectory
rng(45, 'twister');
n = 3;
m = 3;
L = 3;
Tmax = 38; 
%at this seed, infeasible with Tmax < 38
%feasible with Tmax >= 38


LS = lpvsim(n, m, L);

Th_vert = [1 0 0;
           1 0.2 0;
           1 0.3 0.2;
           1.1 0 0;
           1.1 0 0.2;
           1.1 0.1 0.2]';
       
LS.sampler.th = @() sample_th(Th_vert);
traj = LS.sim(Tmax);


%run QMI solver
LP = lpvstab(traj);

out = LP.stab(Th_vert);


%% acquire samples
Nsys = 10;
sys_smp = LS.sample_sys(traj, Nsys);

%% validate samples
if ~out.sol.problem
[valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end