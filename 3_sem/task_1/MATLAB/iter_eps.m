epses = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12; 1e-13; 1e-14; 1e-15];
figure();

%--polyFunc_bisection
iters = [5; 8; 12; 15; 18; 22; 25; 28; 32; 35; 38; 42; 45; 48; 52];
semilogx(epses, iters, 'r', 'LineWidth', 2);
hold on;

%--polyFunc_secant
iters = [4; 5; 5; 6; 6; 7; 7; 7; 7; 8; 8; 8; 8; 8; 8];
semilogx(epses, iters, 'y', 'LineWidth', 2);
hold on;

%--polyFunc_fzero
f = @(x) x.^4 - x - 1;
iters = zeros(size(eps));

for i = 1 : 15
    options = optimset('TolX', epses(i));
    [a, b, c, output] = fzero(f, 0.5, options);
    iters(i) = output.iterations;
end

semilogx(epses, iters, 'g', 'LineWidth', 2);
hold on;

%--transFunc_bisection
iters = [5; 8; 11; 15; 18; 21; 25; 28; 31; 35; 38; 41; 45; 48; 51];
semilogx(epses, iters, 'b', 'LineWidth', 2);
hold on;

%--transFunc_secant
iters = [1; 2; 3; 3; 4; 4; 4; 4; 5; 5; 5; 5; 5; 5; 6];
semilogx(epses, iters, 'm', 'LineWidth', 2);
hold on;

%--transFunc_fzero
f = @(x) x + cos(x);
iters = zeros(size(eps));

for i = 1 : 15
    options = optimset('TolX', epses(i));
    [a, b, c, output] = fzero(f, -1.0, options);
    iters(i) = output.iterations;
end

semilogx(epses, iters, 'c', 'LineWidth', 2);

grid on;
xlabel('epses (accuracy)');
ylabel('number of iterations');
title('Dependence of the number of iterations on accuracy');
legend({'polyFunc bisection', 'polyFunc secant', 'polyFunc fzero', 'transFunc bisection', 'transFunc secant', 'transFunc fzero'}, 'FontSize', 12)
