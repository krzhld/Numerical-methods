clear all;
n = 15; % size of SOLE

% select conditional number
cond_number = 1e5;

% cond number for symmetric matrixes is |lambda_max| / |lambda_min|
lambda = linspace(1, cond_number, n); 
D = diag(lambda);

% get random ortoghonal matrix
[Q, ~] = qr(rand(n));

% ortoghonal matrix don't change eigenvalues, cond(A) = cond(D)
A = Q*D*Q';

suitable_matrix = true;
% checking corner minors
for i = 1 : n
    if det(A(1:i, 1:i)) == 0
        suitable_matrix = false;
    end
end
  
if (suitable_matrix == true)
    % form ground_truth_solution
    x = 10 * rand(n, 1);

    % form column B
    B = A * x;

    fileID = fopen('matrixes.txt', 'w');
    % print conditional number
    fprintf(fileID, '%.5f\n', cond(A));
    % print matrix A
    for i = 1 : n
        for j = 1 : n
            fprintf(fileID, '%.15f ', A(i, j));
        end
        fprintf(fileID, '\n');
    end
    fprintf(fileID, '\n');
    % print column B
    for i = 1 : n
        fprintf(fileID, '%.15f ', B(i));
    end
    fprintf(fileID, '\n\n');
    %print solution
    for i = 1 : n
        fprintf(fileID, '%.15f ', x(i));
    end
    fprintf(fileID, '\n');
    fclose(fileID);
end
    
norm(A, 2)
