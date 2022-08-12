%see if the yalmip robust optimization module can help synthesize the QMIs

LS = lpvsim();

th_bnd = 0.5;
LS.sampler.th = @() (2*rand(LS.L, 1) - 1)*th_bnd;

out = LS.sim(10);

Psi = sample_matrix(out.X, out.U, out.epsilon, out.Th);

delta =1e-5;

%declare variables
P = sdpvar(out.n);
M1 = sdpvar(out.m, out.n, 'full');
M2 = sdpvar(out.m, out.n, 'full');

th = sdpvar(out.L);

