function [Psi] = sample_matrix_lti(X, U, epsilon)
%SAMPLE_MATRIX_LTI Create the data matrix for the QMI under an
%individual-sample noise bound of ||w_t||_2 <= epsilon
%(A, B): [I A' B']' Psi [I A' B'] >= 0
%
%LTI systems only
%
%Inputs:
%   X:          State records
%   U:          Input records
%   epsilon:    L2 noise bound
%
%Outputs:
%   Psi:        Data matrix for the QMI
%get properties of the data
n = size(X, 1);
[m, T] = size(U);

Xp = X(:, 2:end);   %next time step
Xn = X(:, 1:end-1); %current time step

%data matrix
Gamma = [eye(n) Xp; 
         zeros(n), -Xn; 
         zeros(m, n), -U];

%inner matrix with the noise bound
Phi = blkdiag(epsilon^2*T*eye(n), -eye(T));

   
Psi = Gamma*Phi*Gamma';
end

