tabl = zeros(5, 15);
tabl(1, 1:5) = [1.000000e-01 5.452440e-03 7.206218e-01 7 7.304317e-03];
tabl(2, 1:5) = [1.000000e-02 2.237485e-03 2.956838e-01 8 1.228797e-03];
tabl(3, 1:5) = [1.000000e-03 3.770534e-04 4.982539e-02 10 3.486402e-05];
tabl(4, 1:5) = [1.000000e-04 1.548097e-04 2.045705e-02 11 5.876063e-06];
tabl(5, 1:5) = [1.000000e-05 6.356470e-05 8.399615e-03 12 9.905480e-07];
tabl(6, 1:5) = [1.000000e-06 2.610039e-05 3.448970e-03 13 1.669978e-07];
tabl(7, 1:5) = [1.000000e-07 4.400765e-06 5.815270e-04 15 4.747505e-09];
tabl(8, 1:5) = [1.000000e-08 1.807064e-06 2.387894e-04 16 8.006396e-10];
tabl(9, 1:5) = [1.000000e-09 7.420277e-07 9.805314e-05 17 1.351737e-10];
tabl(10, 1:5) = [1.000000e-10 1.251168e-07 1.653320e-05 19 4.035883e-12];
tabl(11, 1:5) = [1.000000e-11 5.137641e-08 6.788988e-06 20 8.526513e-13];
tabl(12, 1:5) = [1.000000e-12 2.109658e-08 2.787747e-06 21 3.410605e-13];
tabl(13, 1:5) = [1.000000e-13 3.557207e-09 4.700568e-07 23 1.989520e-13];
tabl(14, 1:5) = [1.000000e-14 1.011353e-10 1.336425e-08 27 1.705303e-13];
tabl(15, 1:5) = [1.000000e-15 1.011353e-10 1.336425e-08 27 1.705303e-13];

for i = 1 : 15
    epses(i) = tabl(i, 1);
    fact_error(i) = tabl(i, 2);
    discrepancy(i) = tabl(i, 3);
    number_iter(i) = tabl(i, 4);
    lambda_error(i) = tabl(i, 5);
end

figure;
loglog(epses, fact_error, 'r', 'LineWidth', 2);
hold on;
loglog(epses, discrepancy, 'b', 'LineWidth', 2); 
hold on;
loglog(epses(1:10), epses(1:10), '-.b', 'LineWidth', 1);
title('Dependences. Separability = 0.6');
grid on;
xlabel('epses');
ylabel('norm');
legend('Norm of factual error', 'Norm of discrepancy', 'Location', 'southeast', 'FontSize', 14);

figure;

subplot(1, 2, 1);
semilogx(epses, number_iter, 'r', 'LineWidth', 2);
hold on;
title('Dependence number iter on eps. Separability = 0.6');
xlabel('epses');
ylabel('number iterations');
hold on;
grid on;
subplot(1, 2, 2);
loglog(epses, lambda_error, 'b', 'LineWidth', 2);
hold on;
loglog(epses, epses, '-.b', 'LineWidth', 1);
hold on;
title('Dependence lambda error on eps. Separability = 0.6');
xlabel('epses');
ylabel('lambda error');
hold on;
grid on;