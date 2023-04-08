%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('calc1.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol1.txt', 'r');
B = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
subplot(2, 1, 1);
xlabel('x');
ylabel('y');
title('n = 512');
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
fileID = fopen('calc2.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol2.txt', 'r');
B = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

subplot(2, 1, 2);
xlabel('x');
ylabel('y');
title('n = 1024');
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
title('fact error: n = 512');
hold on;

subplot(2, 1, 2);
plot(errors_2_x, errors_2_y, 'LineWidth', 2, 'Color', 'r');
grid on;
xlabel('x');
ylabel('fact error');
title('fact error: n = 1024');
hold on;

figure;
plot(errors_1_x, errors_1_y, 'LineWidth', 2, 'Color', 'r');
hold on;
plot(errors_2_x, errors_2_y, 'LineWidth', 2, 'Color', 'b');
hold on;
legend({'n = 512', 'n = 1024'}, 'Location', 'northeast', 'FontSize', 12);
grid on;
xlabel('x');
ylabel('fact error');
title('fact error');

figure; 
plot(errors_1_x, errors_1_y ./ errors_2_y(1:2:end), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;

%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('dependence_step_error.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
loglog(A(1, 1:11), A(2, 1:11), 'LineWidth', 2, 'Color', 'b');
hold on;
loglog(A(1, 1:11), A(1, 1:11) .^ 3, 'LineWidth', 2, 'Color', 'g');
hold on;
grid on;
xlabel('step');
ylabel('fact error');
title('dependence between step and fact error');
hold on;
legend({'dependence', 'step^3'}, 'Location', 'southeast', 'FontSize', 12);

figure;
semilogx(A(1, 1:11), A(1, 1:11) .^ 3 ./ A(2, 1:11), 'LineWidth', 2, 'Color', 'r');
hold on;
grid on;
xlabel('step');
ylabel('fact error / step^3');
title('dependence between step and fraction');

%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('dependence_pert_error.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
loglog(A(1, 1:12), A(2, 1:12), 'LineWidth', 2, 'Color', 'm');
hold on;
xlabel('perturbation');
ylabel('fact error');
title('dependence between perturbation and fact error (inf norm; step = 0.09)');
grid on;
hold on;
