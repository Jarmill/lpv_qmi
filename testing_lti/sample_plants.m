load('data_test.mat')
rng(30, 'twister');

n = size(Xp, 1);
[m, T] = size(U);

A = sdpvar(n, n, 'full');
B = sdpvar(n, m, 'full');


W = Xp - A*Xn - B*U;
W2 = sum(W.^2, 1);

cons = (W2 <= epsilon^2);

% CA = randn(size(A));
% CB = randn(size(B));

CA = sdpvar(n, n, 'full');
CB = sdpvar(n, m, 'full');

objective = sum(CA.*A, 'all') + sum(CB.*B, 'all');

opts = sdpsettings('solver','mosek');

% sol = optimize(cons, objective, opts);

P = optimizer(cons, objective, opts, [CA, CB], [A,B]);

Ntrial = 20;
sys_rec = cell(Ntrial, 1);
for i = 1:Ntrial
    CA_curr = randn(n, n);
    CB_curr = randn(n, m);
    out_curr = P([CA_curr, CB_curr]);
    sys_rec{i}.A = out_curr(:, 1:n);
    sys_rec{i}.B = out_curr(:, n + (1:m));    
end
% A_rec

% A_rec = value(A)
% B_rec = value(B)