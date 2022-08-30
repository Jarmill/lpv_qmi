x0 = [2; 1.5];
Nsys = length(sys_smp);

LS.epsilon = 0;
% seed = 55;
seed = 200;
% seed = 30;

c = linspecer(2);

sys_aug = [sys_smp;{traj.ground_truth}];
figure(4)
clf
ax1 = subplot(2, 1, 1);
hold on
for i = 1:Nsys+1
    rng(seed, 'twister')

    traj_closed = LS.sim(Tmax, sys_aug{i}, x0);

    if i == Nsys+1
        plot(traj_closed.X(1, :), traj_closed.X(2, :),  'color', c(2, :), 'LineWidth', 3)
        scatter(traj_closed.X(1, :), traj_closed.X(2, :), 10, 'k', 'filled')
    else
        plot(traj_closed.X(1, :), traj_closed.X(2, :), 'color', c(1, :))
        scatter(traj_closed.X(1, :), traj_closed.X(2, :), 20, 'r', 'filled')
    end
    

end

scatter(x0(1), x0(2), 200, 'ok')

xlabel('x_1')
ylabel('x_2')
title('Single Parameter Sequence')

%% multiple parameters
ax2 = subplot(2,1,2);
hold on
seed_num = 12;
for i = 1:Nsys+1
    for j = 1:seed_num
    rng(seed + j, 'twister')

    traj_closed = LS.sim(Tmax, sys_aug{i}, x0);

%     if i == Nsys+1
%         plot(traj_closed.X(1, :), traj_closed.X(2, :), '--', 'color', c(2, :), 'LineWidth', 3)
%     else
        plot(traj_closed.X(1, :), traj_closed.X(2, :), 'color', c(1, :))
        scatter(traj_closed.X(1, :), traj_closed.X(2, :), 10, 'k', 'filled')
%     end
    end

end

title(sprintf('%d Parameter Sequences', seed_num))
scatter(x0(1), x0(2), 200, 'ok')

xlabel('x_1')
ylabel('x_2')


% ploft(