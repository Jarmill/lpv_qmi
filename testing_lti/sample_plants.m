load('data_test.mat')
rng(30, 'twister');

n = size(Xp, 1);
[m, T] = size(U);

A = sdpvar(n, n, 'full');
B = sdpvar(n, m, 'full');


W = Xp - A*Xn - B*U;
W2 = sum(W.^2, 1);

cons = (W2' <= epsilon^2);

CA = randn(size(A));
CB = randn(size(B));

% CA = sdpvar(n, n, 'full');
% CB = sdpvar(n, m, 'full');

objective = sum(CA.*A, 'all') + sum(CB.*B, 'all');

opts = sdpsettings('solver','mosek', 'verbose', 2);

sol = optimize(cons, objective, opts);

% P = optimizer(cons, objective, opts, [CA, CB], [A,B]);

Ntrial = 1;
sys_rec = cell(Ntrial, 1);
i = 1;
% for i = 1:Ntrial
% while i <= Ntrial
%     CA_curr = randn(n, n);
%     CB_curr = randn(n, m);
%     [out_curr, err_curr] = P([CA_curr, CB_curr]);
    
%     if err_curr == 0
%         sys_rec{i}.A = out_curr(:, 1:n);
%         sys_rec{i}.B = out_curr(:, n + (1:m));    
%         sys_rec{i}.W = Xp - sys_rec{i}.A*Xn - sys_rec{i}.B*U;
%         sys_rec{i}.W2 = sum(sys_rec{i}.W.^2, 1);
%         i = i + 1;
%     end
% end


%% see if the recovered matrix satisfies the QMI
A_rec = value(A);
B_rec = value(B);
W_rec = Xp - A_rec*Xn - B_rec*U;
W2_rec = sum(W_rec.^2, 1);

IW = [eye(n); A_rec'; B_rec'];

Psi = sample_matrix_lti([Xn, Xp(:, end)], U, epsilon);

MS = IW' * Psi * IW;
eMS = eig(MS)

% A_rec

