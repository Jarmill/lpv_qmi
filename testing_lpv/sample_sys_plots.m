%run all_2_box_stab_nonsymm.m first
%% simulation setup

x0 = [-2; 1.5];
Nsys = length(sys_smp);

%set up the new sampler
LS2 = lpvsim(n, m, L, 0);
LS2.sampler.th = @() (2*rand(2, 1) - 1) + Th_shift;

%hook up a controller
Pth = get_vertex_interp(Th_vert);
Kth = @(th) K_interp(Pth, out.K, th);
LS2.sampler.u = @(th,x) Kth(th)*x;

% LS2.epsilon = 0;
% seed = 55;
% seed = 200;
seed = 205;
% seed = 30;

c = linspecer(2);

sys_aug = [sys_smp;{traj.ground_truth}];

%% sample trajectories

seed_num = 30;
traj_sim = cell(Nsys+1, seed_num);

for i = 1:(Nsys+1)
    for j = 1:seed_num
        rng(seed + seed_num - 1, 'twister')
        traj_sim{i, j} = LS2.sim(Tmax, sys_aug{i}, x0);
    end
end



%% plot trajectories
figure(4)
clf
ax1 = subplot(2, 1, 1);
hold on
j = 1;
for i = 1:Nsys+1
    rng(seed, 'twister')

%     traj_sim{i, j} = LS2.sim(Tmax, sys_aug{i}, x0);

    if i == Nsys+1
        plot(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :),  'color', c(2, :), 'LineWidth', 3)
        scatter(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 20, c(2, :), 'filled')

%         plot(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :),  'color', c(2, :), 'LineWidth', 3)
%         scatter(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 20, c(2, :), 'filled')
    else
        plot(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 'color', c(1, :))
        scatter(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 10, c(1, :), 'filled')
    end
    

end

scatter(x0(1), x0(2), 200, 'ok')
xlim(ax1, [-2.15, 0.25])
ylim (ax1, [-0.5, 1.6])
xlabel(ax1, '$x_1$','interpreter', 'latex')
ylabel(ax1, '$x_2$', 'interpreter', 'latex')
title(ax1, 'Single Parameter Sequence', 'fontsize', 14)

%% multiple parameters
ax2 = subplot(2,1,2);
hold on
% seed_num = 30;
for i = 1:Nsys+1

%     traj_sim{i, j} = LS2.sim(Tmax, sys_aug{i}, x0);

    for j = 1:seed_num
%     if i == Nsys+1
%         plot(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), '--', 'color', c(2, :), 'LineWidth', 3)
%     else
%         plot(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 'color', c(1, :))
        scatter(traj_sim{i, j}.X(1, :), traj_sim{i, j}.X(2, :), 20, c(1, :), 'filled')
%     end
    end

end

% title(sprintf('%d Parameter Sequences', seed_num))
scatter(x0(1), x0(2), 200, 'ok')

xlim(ax2, [-2.15, 0.25])
ylim (ax2, [-0.5, 1.6])
xlabel(ax2, '$x_1$','interpreter', 'latex')
ylabel(ax2, '$x_2$', 'interpreter', 'latex')
title(ax2, sprintf('%d Parameter Sequences', seed_num), 'fontsize', 14)


% ploft(