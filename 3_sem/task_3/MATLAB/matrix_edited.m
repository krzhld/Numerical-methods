clear all;
clc;

N = 15;
NUMBER_EPSES = 15;

% fixed det
D = diag(linspace(2, 3, N));

writematrix(det(D), 'SOLE.txt', 'Delimiter', 'tab'); 

for iter = 1 : NUMBER_EPSES
    ground_truth = 10 * rand(N, 1);

    w = rand(N, 1);
    Q = eye(N) - 2 * w * w' / (norm(w) * norm(w));
    A = Q' * D * Q;
    
    b = A * ground_truth;
    writematrix(A, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
    writematrix(ground_truth, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
    writematrix(b, 'SOLE.txt', 'Delimiter', 'tab', 'WriteMode', 'append');
 
end


