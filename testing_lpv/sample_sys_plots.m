x0 = [2; 1];
Nsys = length(sys_smp);

figure(4)
clf
hold on
for i = 1:Nsys
    rng(55, 'twister')

    traj_closed = LS.sim(Tmax, sys_smp{i}, x0);

    plot(traj_closed.X(1, :), traj_closed.X(2, :))

end