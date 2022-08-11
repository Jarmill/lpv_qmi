LS = lpvsim();

out = LS.sim(10);

Psi = sample_matrix(out.X, out.U, out.epsilon, out.Th);