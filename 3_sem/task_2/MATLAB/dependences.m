clear all;

fileID = fopen('dependence.txt', 'r');

output = fscanf(fileID, "%e");

for i = 0 : 14
    cond_number(i + 1) = output(3 * i + 1);
    fact_error(i + 1) = output(3 * i + 2);
    discrepancy(i + 1) = output(3 * i + 3);
end

figure;

% subplot(1, 2, 1);
loglog(cond_number, fact_error, 'r', 'LineWidth', 2);
hold on;
grid on;
title('Dependences');
xlabel('conditional number');
ylabel('norm');

% subplot(1, 2, 2);
loglog(cond_number, discrepancy, 'b', 'LineWidth', 2);
hold on;

x = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15];
y = [1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
loglog(x, y, 'g', 'LineWidth', 2);
hold on;

legend('Dependence of factual error on conditional number', 'Dependence of discrepancy on conditional number', 'Bisector', 'Location', 'southeast', 'FontSize', 14)
fclose(fileID);
