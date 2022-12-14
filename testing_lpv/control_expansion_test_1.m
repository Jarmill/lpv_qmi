LS = lpvsim(3, 2, 1);
% LS.L = 1;
out = LS.sim(10);
th_bnd = 1;

%vertices of polytope Theta
Vth = [1, -2];


%% formulate LMI
Y = sdpvar(LS.n, LS.n);
delta = 1e-3;
M1 = sdpvar(LS.m, LS.n);
M2 = sdpvar(LS.m, LS.n);

A = out.ground_truth.A{1};
B = out.ground_truth.B;

S1 = Vth(1)*A*Y + B*M1;
S2 = Vth(2)*A*Y + B*M2;

Y1 = [Y, S1'; S1, Y];
Y2 = [Y, S2'; S2, Y];

cons = [Y>=delta*eye(LS.n); Y1 >= delta*eye(2*LS.n); Y2 >= delta*eye(2*LS.n)];

%% solve LMI and recover controllers
sol = optimize(cons, 0, []);

Y_rec = value(Y);
Y1_rec = value(Y1);
Y2_rec = value(Y2);
M1_rec = value(M1);
M2_rec = value(M2);

K1 = (Y_rec \ M1_rec')';
K2 = (Y_rec \ M2_rec')';

Ac1_rec = Vth(1)*A + B*K1;
Ac2_rec = Vth(2)*A + B*K2;

e1 = eig(Ac1_rec);
e2 = eig(Ac2_rec);

[e1, e2]


%% form and evaluate QMI

IW = [eye(out.n); A'; B'];


V1 = Vth(:, 1);
CT1_bot = [kron(V1*V1', Y_rec), kron(V1, M1_rec)'; kron(V1, M1_rec), M1_rec * (Y_rec \ M1_rec')];

V2 = Vth(:, 2);
CT2_bot = [kron(V2*V2', Y_rec), kron(V2, M2_rec)'; kron(V2, M2_rec), M2_rec * (Y_rec \ M2_rec')];

CT1 = blkdiag(Y_rec, -CT1_bot);
CT2 = blkdiag(Y_rec, -CT2_bot);

Q1 = IW'*CT1*IW;
Q2 = IW'*CT2*IW;




















% Psi = sample_matrix(out.X, out.U, out.epsilon, out.Th);