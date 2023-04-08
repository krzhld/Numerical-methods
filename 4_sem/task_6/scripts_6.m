%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('calc_sol_1.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('exact_sol_1.txt', 'r');
B = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
subplot(2, 1, 1);
xlabel('x');
ylabel('y');
title('step 0.19');
hold on;
n = size(A);
n = n(2);
plot(A(1, 1:n), A(2, 1:n), 'LineWidth', 4, 'Color', 'm');
hold on;
plot(B(1, 1:n), B(2, 1:n), 'LineWidth', 2, 'Color', 'g');
hold on;
grid on;
legend({'calculated', 'solution'}, 'Location', 'southeast', 'FontSize', 12);
errors_1_y = zeros(1, n);
errors_1_y(1, 1:n) = A(2, 1:n) - B(2, 1:n);
errors_1_x = A(1, 1:n);

formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('calc_sol_2.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('exact_sol_2.txt', 'r');
B = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

subplot(2, 1, 2);
xlabel('x');
ylabel('y');
title('step 0.095');
hold on;
n = size(A);
n = n(2);
plot(A(1, 1:n), A(2, 1:n), 'LineWidth', 4, 'Color', 'm');
hold on;
plot(B(1, 1:n), B(2, 1:n), 'LineWidth', 2, 'Color', 'g');
hold on;
grid on;
legend({'calculated', 'solution'}, 'Location', 'southeast', 'FontSize', 12);
errors_2_y = zeros(1, n);
errors_2_y(1, 1:n) = A(2, 1:n) - B(2, 1:n);
errors_2_x = A(1, 1:n);

%----------------------------------------------------
figure;
subplot(2, 1, 1);
plot(errors_1_x, errors_1_y, 'LineWidth', 2, 'Color', 'r');
grid on;
xlabel('x');
ylabel('fact error');
title('fact error: step = 0.19');
hold on;

subplot(2, 1, 2);
plot(errors_2_x, errors_2_y, 'LineWidth', 2, 'Color', 'r');
grid on;
xlabel('x');
ylabel('fact error');
title('fact error: step = 0.095');
hold on;

figure; 
plot(errors_1_x, errors_1_y, 'LineWidth', 2, 'Color', 'r');
hold on;
plot(errors_2_x, errors_2_y, 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;
xlabel('x');
ylabel('fact error');
title('fact error');
legend({'step = 0.19', 'step = 0.095'}, 'Location', 'southeast', 'FontSize', 12);

figure;
plot(errors_1_x(2:1:end), errors_1_y(2:end) ./ errors_2_y(3:2:end), 'LineWidth', 2, 'Color', 'r');
hold on;
grid on;

%----------------------------------------------------
formatSpec = '%f %f %f';
sizeM = [3 Inf];
fileID = fopen('dependence.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
loglog(A(1, 1:18), A(2, 1:18), 'LineWidth', 2, 'Color', 'm');
hold on;
loglog(A(1, 1:18), A(1, 1:18) .^ 3, 'LineWidth', 2, 'Color', 'r');
hold on;
xlabel('step');
ylabel('accuracy');
title('dependence between step and factual accuracy');
hold on;
grid on;
legend({'fact error (inf norm)', 'step^3'}, 'Location', 'southeast', 'FontSize', 12);

figure;
semilogx(A(1, 1:15), A(3, 1:15), 'LineWidth', 2, 'Color', 'r');
hold on;
xlabel('step');
ylabel('constant');
title('dependence between step and constant');
grid on;