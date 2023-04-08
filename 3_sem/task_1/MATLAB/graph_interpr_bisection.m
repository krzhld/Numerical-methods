f = @(x) x.^4 - x - 1;
x = 0.8 : 0.01 : 1.4;
y = f(x);

figure();
grid on;
hold on;
plot(x, y, 'LineWidth', 3);
xlabel('x');
ylabel('y');
title('Graphical interpreation of bisection method');

inter_x = [1.25, 0.875, 1.0625, 1.15625, 1.203125];
inter_y = f(inter_x);

plot(inter_x, inter_y, 'r*');

line([inter_x(1) inter_x(1)], [inter_y(1) 0], 'Color', 'red', 'LineWidth', 2);
line([inter_x(2) inter_x(2)], [inter_y(2) 0], 'Color', 'cyan', 'LineWidth', 2);
line([inter_x(3) inter_x(3)], [inter_y(3) 0], 'Color', 'green', 'LineWidth', 2);
line([inter_x(4) inter_x(4)], [inter_y(4) 0], 'Color', 'blue', 'LineWidth', 2);
line([inter_x(5) inter_x(5)], [inter_y(5) 0], 'Color', 'magenta', 'LineWidth', 2);

legend({'f(x)', 'iteration points', '1 approx', '2 approx', '3 approx', '4 approx', '5 approx'}, 'Location', 'southeast');
