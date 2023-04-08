clear all;
clc;

N = 10;

% fixed eigenvalues
F = linspace(8, 80, N-2);
F(N-1) = 93;
F(N) = 154;

lambda_1 = F(N);
lambda_2 = F(N-1);
lambda_n = F(1);

F = diag(F);
writematrix([lambda_1 lambda_2 lambda_n], 'matrix.txt', 'Delimiter', 'tab'); 

w = rand(N, 1);
Q = eye(N) - 2 * w * w' / (norm(w) * norm(w));
A = Q' * F * Q;

writematrix(A, 'matrix.txt', 'Delimiter', 'tab', 'WriteMode', 'append');

eig_vector_d = linspace(0, 0, N)';
eig_vector_d(N) = 1;
ground_truth = Q' * eig_vector_d;

writematrix(ground_truth, 'matrix.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
