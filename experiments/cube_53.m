%generate the trajectory
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
LP = lpvstab(traj);
% LP.const_K = 1;
out = LP.stab(Th_vert);

out.sol.problem

%H2 control

%% acquire samples

%% validate samples
if ~out.sol.problem
Nsys = 10;
sys_smp = LS.sample_sys(traj, Nsys);

    [valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end

%% plot 
% [k, av] = convhull(Th_vert');
% figure(1)
% Thk = Th_vert(:, k');
% % trisurf(k', Thk(1, :), Thk(2, :), Thk(3, :))
% scatter3(Thk(1, :), Thk(2, :), Thk(3, :), 400, 'filled')
% % trisurf(Th_vert(
%  [A,b,Aeq,beq]=vert2lcon(Th_vert,TOL)