f = @(x) x + cos(x);
x = -2.5 : 0.01 : 1.5;
y = f(x);

figure();
grid on;
hold on;
plot(x, y);
xlabel('x');
ylabel('y');
title('f(x) = x + cos(x)');

fzero(f, -2);
