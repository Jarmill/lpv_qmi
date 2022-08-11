%find the qmi that is satisfied for all plants consistent with data
n = 6;
m = 4;
% T = 15;
T = 40;

%% simulate system
rng(40, 'twister');
sys = drss(n, n, m);

A = sys.a;
B = sys.b;

x0 = randn(n, 1);
epsilon = 0.1;
u_bnd = 1;

X = [x0, zeros(n, T)];
U = zeros(m, T);
W_true = zeros(n, T);

%main simulation loop
xcurr = x0;
for i = 1:T
    %inputs
    ucurr = randn(m, 1)*u_bnd;
    wcurr = randn(n, 1);
    wcurr = epsilon*wcurr/norm(wcurr, 2);
    
    %propagation
    xnext = A*xcurr + B*ucurr + wcurr;
    
    %storage
    X(:, i+1) = xnext;
    U(:, i) = ucurr;
    W_true(:, i) = wcurr;
    xcurr = xnext;
end

%compute bounds
Xp = X(:, 2:end);
Xn = X(:, 1:end-1);

%recovered noise values
W = Xp - A*Xn - B*U;

%% set up QMI constraints
IW = [eye(n); W'];
Phi = blkdiag(epsilon^2*T*eye(n), -eye(T));


X = IW'*Phi*IW;
eX = eig(X);

Gamma = [eye(n) Xp; zeros(n), -Xn; zeros(m, n), -U];

Psi = Gamma*Phi*Gamma';


%all plants consistent with data must have Mtrue >= 0

Strue = [eye(n); A'; B'];
Mtrue = Strue'*Psi*Strue


% sysfull = ss(A, B, eye(n), zeros(n, m), 1);