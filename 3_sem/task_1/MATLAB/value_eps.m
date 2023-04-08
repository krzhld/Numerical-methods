epses = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12; 1e-13; 1e-14; 1e-15];
figure();

%--polyFunc_bisection
polyFunc = [1.180095e-01; 2.584824e-02; 4.325295e-04; 2.122154e-04; 1.945875e-05; 5.282001e-07; 1.171568e-07; 1.729253e-08; 1.194251e-09; 1.438663e-10; 1.256772e-11; 1.080913e-12; 5.506706e-14; 2.264855e-14; 2.220446e-15];
loglog(epses, polyFunc, 'r', 'LineWidth', 2);
hold on;

%--polyFunc_secant
polyFunc = [1.263697e-02; 2.967623e-04; 2.967623e-04; 8.471037e-07; 8.471037e-07; 5.704770e-11; 5.704770e-11; 5.704770e-11; 5.704770e-11; 4.440892e-16; 4.440892e-16; 4.440892e-16; 4.440892e-16; 4.440892e-16; 4.440892e-16];
loglog(epses, polyFunc, 'y', 'LineWidth', 2);
hold on;

%--polyFunc_fzero
f = @(x) x.^4 - x - 1;
polyFunc = zeros(size(eps));

for i = 1 : 15
    options = optimset('TolX', epses(i));
    polyFunc(i) = abs(f(fzero(f, 0.5, options)));
end

loglog(epses, polyFunc, 'g', 'LineWidth', 2);
hold on;

%--transFunc_bisection
transFunc = [3.387937e-02; 1.345150e-03; 5.281584e-04; 3.362535e-05; 1.703583e-06; 6.905383e-07; 4.213050e-08; 1.512334e-09; 7.329977e-10; 2.369882e-12; 3.718692e-12; 8.659740e-14; 3.907985e-14; 2.553513e-15; 3.330669e-16];
loglog(epses, transFunc, 'b', 'LineWidth', 2);

%--transFunc_secant
transFunc = [4.012602e-02; 1.974111e-03; 1.025799e-05; 1.025799e-05; 2.669348e-09; 2.669348e-09; 2.669348e-09; 2.669348e-09; 3.552714e-15; 3.552714e-15; 3.552714e-15; 3.552714e-15; 3.552714e-15; 3.552714e-15; 0.000000e+00];
for i = 1 : 15
    if transFunc(i) == 0.000000e+00
        transFunc(i) = 1.000000e-17;
    end
end
    
loglog(epses, transFunc, 'm', 'LineWidth', 2);

%--transFunc_fzero
f = @(x) x + cos(x);
transFunc = zeros(size(eps));

for i = 1 : 15
    options = optimset('TolX', epses(i));
    transFunc(i) = abs(f(fzero(f, -1.0, options)));
end

loglog(epses, transFunc, 'c', 'LineWidth', 2);

ylim([1.000000e-17 1.000000e1])
grid on;
xlabel('epses (accuracy)');
ylabel('module of value of function near root');
title('Dependence of module of value near root on accuracy');
legend({'polyFunc bisection', 'polyFunc secant', 'polyFunc fzero', 'transFunc bisection', 'transFunc secant', 'transFunc fzero'}, 'FontSize', 12, 'Location', 'SouthEast')
