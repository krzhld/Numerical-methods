#include <stdio.h>
#include <math.h>

double PolyFunc(double x) {
	double res = pow(x, 4.0) - x - 1;
	return res;
}

double TransFunc(double x) {
	double res = x + cos(x);
	return res;
}

double BisectionMethod(double(*func)(double), int* numIter, double a, double b, double eps) {
	int iter = 0;

	double func_a = func(a);
	double c;
	while (b - a >= eps) {
		c = (a + b) / 2;
		if (func(c) * func_a < 0)
			b = c;
		else
			a = c;
		++iter;
	}

	if (numIter != NULL)
		*numIter = iter;

	return ((b + a) / 2);
}

/* m1 - min f', M2 - max f'' */
double SecantMethod(double(*func)(double), int* numIter, double x1, double x2, double m1, double M2, double eps) {
	int iter = 0;

	double beta = (sqrt(5.0) + 1) / 2.0;
	double C = pow((2 * m1 / M2), 1 / beta);

	double denominator;
	double x3;
	do {
		denominator = func(x2) - func(x1);
		x3 = -func(x2) * (x2 - x1) / denominator + x2;
		x1 = x2;
		x2 = x3;
		++iter;
	} while (pow(fabs(x2 - x1), beta) >= eps * C);

	if (numIter != NULL)
		*numIter = iter;

	return x3;
}

int main() {
	double eps = 1e-15;
	
	// Деление пополам. Корни.
	printf("Bisection method\n");
	//printf("Polynominal function, negative root: %lf\n", BisectionMethod(PolyFunc, NULL, -2, -0.5, eps));
	printf("Polynominal function, positive root: %lf\n", BisectionMethod(PolyFunc, NULL, 0.7, 3, eps));
	printf("Transcendental function, root: %lf\n", BisectionMethod(TransFunc, NULL, -1.55, -0.1, eps));

	// Метод секущих. Корни. Начальные приближения взяли по усл Фурье
	printf("\nSecant method\n");
	//printf("Polynominal function, negative root:%lf\n", SecantMethod(PolyFunc, NULL, -1.5, -1, eps));
	/* PolyFunc(x) = x^4 - x - 1, f' = 4x^3 - 1, f'' = 12x^2, [a; b] = [0.7; 3], m1 = 0.372, M2 = 108 */
	double m1_polyFunc = 0.372;
	double M2_polyFunc = 108.0;
	printf("Polynominal function, positive root: %lf\n", SecantMethod(PolyFunc, NULL, 1.5, 2, m1_polyFunc, M2_polyFunc, eps));
	/* TransFunc(x) = x + cos(x), f' = 1 - sin(x), f'' = -cos(x), [a; b] = [-1.55; -0.1], m1 = 1, M2 = 1 */
	double m1_transFunc = 1;
	double M2_transFunc = 1;
	printf("Transcendental function, root: %lf\n", SecantMethod(TransFunc, NULL, -1.5, -1, m1_transFunc, M2_transFunc, eps));
	
	double epses[15] = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15 };
	
	double polyFunc_bisection[15];
	double transFunc_bisection[15];
	int iter_polyFunc_bisection[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int iter_transFunc_bisection[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	double polyFunc_secant[15];
	double transFunc_secant[15];
	int iter_polyFunc_secant[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int iter_transFunc_secant[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	double curRoot_polyFunc, curRoot_transFunc;

	for (int i = 0; i < 15; i++) {
		curRoot_polyFunc = BisectionMethod(PolyFunc, iter_polyFunc_bisection + i, 0.7, 3, epses[i]);
		curRoot_transFunc = BisectionMethod(TransFunc, iter_transFunc_bisection + i, -2.0, 0.0, epses[i]);
		polyFunc_bisection[i] = fabs(PolyFunc(curRoot_polyFunc));
		transFunc_bisection[i] = fabs(TransFunc(curRoot_transFunc));

		curRoot_polyFunc = SecantMethod(PolyFunc, iter_polyFunc_secant + i, 1.5, 2, m1_polyFunc, M2_polyFunc, epses[i]);
		curRoot_transFunc = SecantMethod(TransFunc, iter_transFunc_secant + i, -1.5, -1, m1_transFunc, M2_transFunc, epses[i]);
		polyFunc_secant[i] = fabs(PolyFunc(curRoot_polyFunc));
		transFunc_secant[i] = fabs(TransFunc(curRoot_transFunc));
	}

	printf("\n\nAbsolute error\n\n");

	printf("Polynominal function, bisection method\n");
	for (int i = 0; i < 15; i++) {
		printf("%le; ", polyFunc_bisection[i]);
	}

	printf("\n\nTranscendental function, bisection method\n");
	for (int i = 0; i < 15; i++) {
		printf("%le; ", transFunc_bisection[i]);
	}

	printf("\n\nPolynominal function, secant method\n");
	for (int i = 0; i < 15; i++) {
		printf("%le; ", polyFunc_secant[i]);
	}

	printf("\n\nTranscendental function, secant method\n");
	for (int i = 0; i < 15; i++) {
		printf("%le; ", transFunc_secant[i]);
	}
	
	
	// Number of iterations for finding root with bisection method
	
	printf("\n\n\nNumber of iterations\n\n");

	printf("Polynominal function, bisection method\n");
	for (int i = 0; i < 15; i++) {
		printf("%d; ", iter_polyFunc_bisection[i]);
	}

	printf("\n\nTranscendental function, bisection method\n");
	for (int i = 0; i < 15; i++) {
		printf("%d; ", iter_transFunc_bisection[i]);
	}

	printf("\n\nPolynominal function, secant method\n");
	for (int i = 0; i < 15; i++) {
		printf("%d; ", iter_polyFunc_secant[i]);
	}

	printf("\n\nTranscendental function, secant method\n");
	for (int i = 0; i < 15; i++) {
		printf("%d; ", iter_transFunc_secant[i]);
	}
	printf("\n");
	
	return 0;
}