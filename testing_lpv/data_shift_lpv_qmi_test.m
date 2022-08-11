%find the qmi that is satisfied for all plants consistent with data
n = 4;
m = 3;
T = 10;
L = 2;

%% simulate system
rng(40, 'twister');
A = cell(L, 1);
for i = 1:L
    sys = drss(n, n, m);
    A{i} = sys.a;
    if i == 1
        B = sys.b;
    end
end


x0 = randn(n, 1);
epsilon = 0.1;
u_bnd = 1;
th_bnd = 1;

X = [x0, zeros(n, T)];
U = zeros(m, T);
W_true = zeros(n, T);
Th = zeros(L, T);
%main simulation loop
xcurr = x0;
for i = 1:T
    %inputs
    ucurr = randn(m, 1)*u_bnd;
    wcurr = randn(n, 1);
    wcurr = epsilon*wcurr/norm(wcurr, 2);
    thcurr = (2*rand(L, 1) - 1)*th_bnd;
    
    %propagation
    xnext = B*ucurr + wcurr;
    for k = 1:L
        xnext = xnext + A{k}*xcurr*thcurr(k);
    end
    %     xnext = A*xcurr + B*ucurr + wcurr;
    
    %storage
    X(:, i+1) = xnext;
    U(:, i) = ucurr;
    W_true(:, i) = wcurr;
    Th(:, i) = thcurr;
    xcurr = xnext;
end

%compute bounds
Xp = X(:, 2:end);
Xn = X(:, 1:end-1);

%recovered noise values
W = Xp - B*U;
for k = 1:L
    W = W - (A{k}*Xn).*Th(k, :);
end

W2 = sum(W.^2, 1)
%% set up QMI constraints
IW = [eye(n); W'];
Phi = blkdiag(epsilon^2*T*eye(n), -eye(T));
% 
% 
MS = IW'*Phi*IW;
eMS = eig(MS);
% 
Gamma = zeros(n+T, n*(L+1) + m);
Gamma(1:n, n+1:end) = Xp;
for k = 1:L
    curr_ind = (1:n) + n*k;
    Gamma(curr_ind, n+1:end) = -Xn.*Th(k, :);
end
Gamma(end-m:end, n+1:end) = -U;
% Gamma = [eye(n) Xp; zeros(n), -Xn; zeros(m, n), -U];
% 
Psi = Gamma*Phi*Gamma';
% 
% 
% %all plants consistent with data must have Mtrue >= 0
% 
% Strue = [eye(n); A'; B'];
% Mtrue = Strue'*Psi*Strue;
% 
% eig(Mtrue)
% 
% sysfull = ss(A, B, eye(n), zeros(n, m), 1);