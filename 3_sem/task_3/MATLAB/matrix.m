clear all;
clc;

N = 15;
NUMBER_EPSES = 15;

% fixed det
F = diag(linspace(10, 11, N));

writematrix(det(F), 'SOLE.txt', 'Delimiter', 'tab'); 

for iter = 1 : NUMBER_EPSES
    ground_truth = 10 * rand(N, 1);
    
    for i = 1 : N
        F(i, i+1 : N) = 0.01 * rand(1, N - i);
    end

    w = rand(N, 1);
    Q = eye(N) - 2 * w * w' / (norm(w) * norm(w));
    % ||A|| <= ||P'||*||B||*||P|| = ||B||
    A = Q' * F * Q;
    % norm(eye(N) - A / norm(A, Inf), Inf);
    
    b = A * ground_truth;
    writematrix(A, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
    writematrix(ground_truth, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
    writematrix(b, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
 
end


