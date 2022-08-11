%try to synthesize a stabilizing controller

load('data_test.mat')
rng(30, 'twister');

%set up properties
n = size(X, 1);
[m, T] = size(U);

Psi = sample_matrix(X, U, epsilon);
delta = 1e-6;


%% declare variables
alpha = sdpvar(1,1);
beta = sdpvar(1,1);

P = sdpvar(n, n);
L = sdpvar(m, n);


%set up QMI
CTRL_TOP = P - beta*eye(n);
CTRL_BOT = [-P, -L', zeros(n, n);
            -L, zeros(m, m), L;
            zeros(n, n), L', P];
CTRL = blkdiag(CTRL_TOP, CTRL_BOT);

DATA = blkdiag(Psi, zeros(n,n));


cons = [alpha>=0;
        beta>=delta;
        P >= delta*eye(n);
        CTRL - alpha*DATA >= 0];        
    

%% solve program
opts = sdpsettings('solver', 'mosek');
sol = optimize(cons, 0, opts);


%% recover
P_rec = value(P);
L_rec = value(L);
K_rec = (P_rec \ L_rec')';

%A + BK should be stable 
