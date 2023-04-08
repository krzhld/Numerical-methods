fileID = fopen('dependence_edited.txt', 'r');

output = fscanf(fileID, "%e");
fclose(fileID);

for i = 0 : 14
    epses(i + 1) = output(4 * i + 1);
    fact_error(i + 1) = output(4 * i + 2);
    discrepancy(i + 1) = output(4 * i + 3);
    number_iter(i + 1) = output(4 * i + 4);
end

figure;

subplot(1, 2, 1);
loglog(epses, fact_error, 'r', 'LineWidth', 2);
hold on;
grid on;
loglog(epses, discrepancy, 'b', 'LineWidth', 2);
title('Dependences');
xlabel('eps');
ylabel('norm');
legend('Dependence of factual error on accuracy', 'Dependence of discrepancy on accuracy', 'Location', 'southeast', 'FontSize', 14)

subplot(1, 2, 2);
semilogx(epses, number_iter, 'g', 'LineWidth', 2);
title('Dependence');
xlabel('eps');
ylabel('number iterations');
hold on;
grid on;
legend('Dependence of number iterations on accuracy', 'Location', 'southeast', 'FontSize', 14);
