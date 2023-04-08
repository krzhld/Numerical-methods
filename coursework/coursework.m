%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('sol_exact_1.txt', 'r');
exact_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol_rungekutta_1.txt', 'r');
rungekutta_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol_predcorr_1.txt', 'r');
predcorr_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
xlabel('x');
ylabel('y');
title('step 0.11875 (n = 16)');
hold on;
plot(exact_sol(1, 1:end), exact_sol(2, 1:end), 'LineWidth', 4, 'Color', 'r');
hold on;
plot(rungekutta_sol(1, 1:end), rungekutta_sol(2, 1:end), 'LineWidth', 3, 'Color', 'g');
hold on;
plot(predcorr_sol(1, 1:end), predcorr_sol(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;
legend({'solution', 'runge kutta', 'predcorr'}, 'Location', 'southeast', 'FontSize', 12);

formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('sol_exact_2.txt', 'r');
exact_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol_rungekutta_2.txt', 'r');
rungekutta_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('sol_predcorr_2.txt', 'r');
predcorr_sol = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
xlabel('x');
ylabel('y');
title('step 0.059375 (n = 32)');
hold on;
plot(exact_sol(1, 1:end), exact_sol(2, 1:end), 'LineWidth', 4, 'Color', 'r');
hold on;
plot(rungekutta_sol(1, 1:end), rungekutta_sol(2, 1:end), 'LineWidth', 3, 'Color', 'g');
hold on;
plot(predcorr_sol(1, 1:end), predcorr_sol(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;
legend({'solution', 'runge kutta', 'predcorr'}, 'Location', 'southeast', 'FontSize', 12);


%----------------------------------------------------
formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('errors_rungekutta_1.txt', 'r');
rungekutta_err_1 = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('errors_predcorr_1.txt', 'r');
predcorr_err_1 = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
xlabel('x');
ylabel('fact error');
title('fact error: step 0.11875 (n = 16)');
hold on;
plot(rungekutta_err_1(1, 1:end), rungekutta_err_1(2, 1:end), 'LineWidth', 2, 'Color', 'r');
hold on;
plot(predcorr_err_1(1, 1:end), predcorr_err_1(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
legend({'runge kutta', 'predcorr'}, 'Location', 'southeast', 'FontSize', 12);
grid on;


formatSpec = '%f %f';
sizeM = [2 Inf];
fileID = fopen('errors_rungekutta_2.txt', 'r');
rungekutta_err_2 = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);
fileID = fopen('errors_predcorr_2.txt', 'r');
predcorr_err_2 = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
xlabel('x');
ylabel('fact error');
title('fact error: step 0.059375 (n = 32)');
hold on;
plot(rungekutta_err_2(1, 1:end), rungekutta_err_2(2, 1:end), 'LineWidth', 2, 'Color', 'r');
hold on;
plot(predcorr_err_2(1, 1:end), predcorr_err_2(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
legend({'runge kutta', 'predcorr'}, 'Location', 'southeast', 'FontSize', 12);
grid on;


figure;
plot(rungekutta_err_1(1, 1:end), rungekutta_err_1(2, 1:end), 'LineWidth', 2, 'Color', 'r');
hold on;
plot(predcorr_err_1(1, 1:end), predcorr_err_1(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
plot(rungekutta_err_2(1, 1:end), rungekutta_err_2(2, 1:end), 'LineWidth', 2, 'Color', 'r');
hold on;
plot(predcorr_err_2(1, 1:end), predcorr_err_2(2, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;
xlabel('x');
ylabel('fact error');
legend({'runge kutta 1', 'predcorr 1', 'runge kutta 2', 'predcorr 2'});

% figure;
% plot(rungekutta_err_1(1, 1:end), predcorr_err_1(2, 1:end) ./ rungekutta_err_1(2, 1:end), 'LineWidth', 2, 'Color', 'r');
% hold on;
% grid on;
% 
% figure;
% plot(rungekutta_err_2(1, 1:end), predcorr_err_2(2, 1:end) ./ rungekutta_err_2(2, 1:end), 'LineWidth', 2, 'Color', 'r');
% hold on;
% grid on;

figure;
plot(rungekutta_err_1(1, 1:end), rungekutta_err_1(2, 1:end) ./ rungekutta_err_2(2, 1:2:end), 'LineWidth', 2, 'Color', 'r');
hold on;
plot(predcorr_err_1(1, 1:end), predcorr_err_1(2, 1:end) ./ predcorr_err_2(2, 1:2:end), 'LineWidth', 2, 'Color', 'b');
hold on;
grid on;

%----------------------------------------------------
formatSpec = '%d %f %f %f %f';
sizeM = [5 Inf];
fileID = fopen('dependences.txt', 'r');
A = fscanf(fileID, formatSpec, sizeM);
fclose(fileID);

figure;
loglog(A(1, 1:end), A(2, 1:end), 'LineWidth', 2, 'Color', 'r');
hold on;
loglog(A(1, 1:end), A(3, 1:end), 'LineWidth', 2, 'Color', 'b');
hold on;
legend({'runge kutta', 'predcorr'}, 'Location', 'southeast', 'FontSize', 12);
hold on;
grid on;
xlabel('number of split segments');
ylabel('norm of fact error');
title('inf norm of fact error');
hold on;
