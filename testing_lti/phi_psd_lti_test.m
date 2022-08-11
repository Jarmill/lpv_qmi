rng(20, 'twister')

%individual sample noise bounds
%case (ii) of page 3 of https://arxiv.org/pdf/2203.12959.pdf

d = 4;
T = 8;

epsilon = 2;

W = randn(d, T);

W = epsilon*normalize(W, 1, 'norm', 2);


IW = [eye(d); W'];
Phi = blkdiag(epsilon^2*T*eye(d), -eye(T));
% 
X = IW'*Phi*IW;
eX = eig(X)

% M = [epsilon^2*eye(N), W; W', eye(d)]