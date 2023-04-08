%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('step1_calculated_solution.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('step1_exact_solution.txt', 'r');
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
fileID = fopen('step2_calculated_solution.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('step2_exact_solution.txt', 'r');
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
plot(errors_1_x(2:1:end), errors_1_y(2:end) ./ errors_2_y(3:2:end), 'LineWidth', 2, 'Color', 'r');
hold on;
grid on;


%----------------------------------------------------
formatSpec = '%f %f %d';
sizeM = [3 Inf];
fileID = fopen('dependences.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
subplot(1, 2, 1);
loglog(A(1, 1:14), A(2, 1:14), 'LineWidth', 2, 'Color', 'm');
hold on;
loglog(A(1, 1:14), A(1, 1:14), 'LineWidth', 2, 'Color', 'r');
hold on;
xlabel('eps');
ylabel('error');
title('dependence between eps and error');
hold on;
grid on;
legend({'fact error (inf norm)', 'bisector'}, 'Location', 'southeast', 'FontSize', 12);

subplot(1, 2, 2);
semilogx(A(1, 1:14), A(3, 1:14), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;
xlabel('eps');
ylabel('number of iterations');
title('dependence between eps and number of iterations');
hold on;

%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('perturbation.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
loglog(A(1, 1:12), A(2, 1:12), 'LineWidth', 2, 'Color', 'm');
hold on;
xlabel('perturbation');
ylabel('fact error');
title('dependence between perturbation and fact error (inf norm; eps = 1e-4)');
grid on;
hold on;
