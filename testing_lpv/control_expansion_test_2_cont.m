%continuous-time indexing
LS = lpvsim(4, 3, 2);
out = LS.sim(10);
th_bnd = 1;

%vertices of polytope Theta
Vth = [1 0 1;
       0 1 1];
   
   n = LS.n;
   m = LS.m;
   L = LS.L;


%% formulate LMI
Y = sdpvar(LS.n, LS.n);
delta = 1e-3;
M1 = sdpvar(LS.m, LS.n);
M2 = sdpvar(LS.m, LS.n);
M3 = sdpvar(LS.m, LS.n);

A1 = out.ground_truth.A{1};
A2 = out.ground_truth.A{2};
B = out.ground_truth.B;

S1 = (Vth(1, 1)*A1 + Vth(2, 1)*A2)*Y+ B*M1;
S2 = (Vth(1, 2)*A1 + Vth(2, 2)*A2)*Y+ B*M2;
S3 = (Vth(1, 3)*A1 + Vth(2, 3)*A2)*Y+ B*M3;

Y1 = -(S1' + S1)/2;
Y2 = -(S2' + S2)/2;
Y3 = -(S3' + S3)/2;

%P_rec = inv(Y_rec)

cons = [Y>=delta*eye(LS.n); Y1 >= delta*eye(LS.n); Y2 >= delta*eye(LS.n)...
    Y3 >= delta*eye(LS.n)];

%% solve LMI and recover controllers
sol = optimize(cons, 0, sdpsettings('solver', 'mosek'));

Y_rec = value(Y);
Y1_rec = value(Y1);
Y2_rec = value(Y2);
Y3_rec = value(Y3);
M1_rec = value(M1);
M2_rec = value(M2);
M3_rec = value(M3);

K1 = (Y_rec \ M1_rec')';
K2 = (Y_rec \ M2_rec')';
K3 = (Y_rec \ M3_rec')';

% Ac1_rec = Vth(1)*A + B*K1;
% Ac2_rec = Vth(2)*A + B*K2;
Ac1_rec = (Vth(1, 1)*A1 + Vth(2, 1)*A2)+ B*K1;
Ac2_rec = (Vth(1, 2)*A1 + Vth(2, 2)*A2)+ B*K2;
Ac3_rec = (Vth(1, 3)*A1 + Vth(2, 3)*A2)+ B*K3;



e1 = eig(Ac1_rec);
e2 = eig(Ac2_rec);
e3 = eig(Ac3_rec);

[e1, e2, e3]


%% form and evaluate QMI

CT1_row = [Vth(1, 1)*Y_rec, Vth(2, 1)*Y_rec, M1_rec'];
CT1 = -[zeros(n), CT1_row; CT1_row', zeros(L*n + m)];


CT2_row = [Vth(1, 2)*Y_rec, Vth(2, 2)*Y_rec, M2_rec'];
CT2 = -[zeros(n), CT2_row; CT2_row', zeros(L*n + m)];

CT3_row = [Vth(1, 3)*Y_rec, Vth(2, 3)*Y_rec, M3_rec'];
CT3 = -[zeros(n), CT3_row; CT3_row', zeros(L*n + m)];


% CT = [


IW = [eye(out.n); A1'; A2'; B'];

% 
% V1 = Vth(:, 1);
% CT1_bot = [kron(V1*V1', Y_rec), kron(V1, M1_rec'); kron(V1', M1_rec), M1_rec * (Y_rec \ M1_rec')];
% 
% V2 = Vth(:, 2);
% CT2_bot = [kron(V2*V2', Y_rec), kron(V2, M2_rec'); kron(V2', M2_rec), M2_rec * (Y_rec \ M2_rec')];
% 
% V3 = Vth(:, 3);
% CT3_bot = [kron(V3*V3', Y_rec), kron(V3, M3_rec'); kron(V3', M3_rec), M3_rec * (Y_rec \ M3_rec')];
% 
% 
% CT1 = blkdiag(Y_rec, -CT1_bot);
% CT2 = blkdiag(Y_rec, -CT2_bot);
% CT3 = blkdiag(Y_rec, -CT3_bot);
% 
Q1 = IW'*CT1*IW;
Q2 = IW'*CT2*IW;
Q3 = IW'*CT3*IW;





%% test indexing


% CT1_schur = [kron(V1*V1', Y_rec), kron(V1, M1_rec'), zeros(n*L, n);
%     kron(V1', M1_rec), zeros(m,m), M1_rec;
%     zeros(n, L*n), M1_rec', Y_rec];
% 











% Psi = sample_matrix(out.X, out.U, out.epsilon, out.Th);