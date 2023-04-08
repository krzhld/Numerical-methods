shifts = [0 5 10 15 20 25 30 35 40 45 50 55 60];
iters = [21 20 19 18 18 17 16 15 14 14 13 16 19];

plot(shifts, iters, '-bo', 'LineWidth', 2);
grid on;
title('Dependence number iterations on shift. EPS = 10^{-6}. Separability = 0.6');
xlabel('shift (optimal shift 50.5)');
ylabel('number of iterations');
legend('scalar multiplication method with shift', 'FontSize', 10);
