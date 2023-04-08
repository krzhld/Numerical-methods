f = @(x) x.^4 - x - 1;
x = -2.5 : 0.01 : 2.5;
y = f(x);

figure();
grid on;
hold on;
plot(x, y);
xlabel('x');
ylabel('y');
title('f(x) = x^4 - x - 1');

fzero(f, -1);
fzero(f, 1);
