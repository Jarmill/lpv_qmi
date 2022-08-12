sdpvar x w
F = [x+w <= 1];
W = [-0.5 <= w <= 0.5, uncertain(w)];
objective = -x;
sol = optimize(F + W,objective);