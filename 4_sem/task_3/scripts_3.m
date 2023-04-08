n = 12;


fileID = fopen('error_and_iter_on_eps.txt', 'r');
formatSpec = '%f %f %d %f';
sizeA = [4 n];
A = fscanf(fileID, formatSpec, sizeA);
fclose(fileID);


figure;
loglog(A(1, 1:n), A(2, 1:n), 'LineWidth', 3, 'Color', 'b');
grid on;
hold on;
loglog(A(1, 1:n), A(1, 1:n), 'LineWidth', 3, 'Color', 'r');
hold on;
loglog(A(1, 1:n), A(2, 1:n), 'mo', 'MarkerSize', 7, 'LineWidth', 3);
title('Fact error(eps)');
xlabel('Accuracy');
ylabel('Error');
legend({'Fact error', 'Bisector'}, 'Location', 'southeast', 'Fontsize', 12);


figure;
semilogx(A(1, 1:n), A(3, 1:n), 'LineWidth', 3);
grid on;
hold on;
semilogx(A(1, 1:n), A(3, 1:n), 'mo', 'MarkerSize', 7, 'LineWidth', 3);
title('Number iter(eps)');
xlabel('Accuracy');
ylabel('Number iter (degree of 2)');


figure;
loglog(A(4, 1:n), A(2, 1:n), 'LineWidth', 3);
grid on;
hold on;
loglog(A(4, 1:n), A(2, 1:n), 'mo', 'MarkerSize', 7, 'LineWidth', 3);
title('Fact error(rank)');
xlabel('Rank of partition');
ylabel('Error');
hold on;

for i = 1:11
constant = (A(2, i + 1) - A(2, i)) / (A(4, i + 1) .^ 4 - A(4, i) .^ 4)
end
