#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

inline double Solution(double x) {
	return exp(2 * x);
}

inline double DerivativeSolution(double x) {
	return 2 * exp(2 * x);
}

inline double f1(double x, double y1, double y2) {
	return y2;
}

double f2_inhomogeneous(double x, double y1, double y2) {
	/*if (x == 0)
		return 4;*/
	double res = exp(2 * x) - 3 * y1 + (2 * x + 1) * y2;
	res /= x;
	return res;
}

double f2_homogeneous(double x, double y1, double y2) {
	/*if (x == 0)
		return 0;*/
	double res = - 3 * y1 + (2 * x + 1) * y2;
	res /= x;
	return res;
}

double GetFactError(double a, double b, int n, vector<double>& grid_function, double(*sol)(double)) {
	double h = (b - a) / n;
	double x = a;
	double result = 0;
	double temp;
	for (int i = 0; i <= n; i++) {
		temp = abs(grid_function[i] - sol(x));
		if (temp > result)
			result = temp;
		x += h;
	}
	return result;
}

// method has 3-rd order coef 1/2
vector<double> RungeKuttaMethod(double a, double b, int n, double y0_a, double y1_a, double(*func1)(double, double, double), double(*func2)(double, double, double)) {
	vector<double> grid_function;
	double h = (b - a) / n;
	double k1, k2, k3;
	double q1, q2, q3;

	double x = a;
	double y0 = y0_a;
	double y1 = y1_a;
	for (int i = 0; i <= n; i++) {
		grid_function.push_back(y0);

		k1 = func1(x, y0, y1);
		q1 = func2(x, y0, y1);

		k2 = func1(x + h / 2, y0 + h * k1 / 2, y1 + h * q1 / 2);
		q2 = func2(x + h / 2, y0 + h * k1 / 2, y1 + h * q1 / 2);

		k3 = func1(x + h, y0 - h * k1 + 2 * h * k2, y1 - h * q1 + 2 * h * q2);
		q3 = func2(x + h, y0 - h * k1 + 2 * h * k2, y1 - h * q1 + 2 * h * q2);

		y0 = y0 + h * (k1 + 4 * k2 + k3) / 6;
		y1 = y1 + h * (q1 + 4 * q2 + q3) / 6;

		x += h;
	}

	return grid_function;
}

vector<double> SolveBoundaryProblem(double a, double b, int n, double alpha_0, double alpha_1, double A, double B) {
	// solution of inhomogeneous equation
	vector<double> u = RungeKuttaMethod(a, b, n, 0, 0, f1, f2_inhomogeneous);
	// linear independent solutions of homogeneous equation
	vector<double> v = RungeKuttaMethod(a, b, n, 1, 0, f1, f2_homogeneous);
	vector<double> w = RungeKuttaMethod(a, b, n, 0, 1, f1, f2_homogeneous);

	// solve 2*2 SOLE 
	double determinant = alpha_0 * w[n] - alpha_1 * v[n];
	if (determinant == 0)
		return { -1.0 };

	double c1 = (A * w[n] - alpha_1 * (B - u[n])) / determinant;
	double c2 = (alpha_0 * (B - u[n]) - A * v[n]) / determinant;

	vector<double> solution(n + 1);
	for (int i = 0; i <= n; i++)
		solution[i] = u[i] + c1 * v[i] + c2 * w[i];

	return solution;
}

void WriteGridFunctionToFile(double a, double b, int n, string file_name, vector<double>& grid_function) {
	double h = (b - a) / n;
	double x = a;

	ofstream outfile(file_name);
	outfile.precision(15);
	for (int i = 0; i <= n; i++) {
		x = a + i * h;
		outfile << x << " " << grid_function[i] << endl;
	}
	outfile.close();
}

vector<double> SolutionToGridFunction(double(*f)(double), double a, double b, int n) {
	vector<double> grid_function;
	double h = (b - a) / n;

	for (int i = 0; i <= n; i++)
		grid_function.push_back(f(a + i * h));

	return grid_function;
}

void FormDependenceStepError(double a, double b, double(*f)(double), double alpha_0, double alpha_1, double A, double B) {
	double h;
	vector<double> solution;

	ofstream outfile("dependence_step_error.txt");
	outfile.precision(15);

	for (int n = 2; n < 2049; n *= 2) {
		solution = SolveBoundaryProblem(a, b, n, alpha_0, alpha_1, A, B);
		h = (b - a) / n;

		outfile << h << " " << GetFactError(a, b, n, solution, f) << endl;
	}
	outfile.close();
}

void FormDependencePerturbationError(double a, double b, int n, double(*f)(double), double alpha_0, double alpha_1, double A, double B) {
	auto perturbation = { 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9 };
	double pert_A;
	double pert_B;
	vector<double> solution;

	ofstream outfile("dependence_pert_error.txt");
	outfile.precision(15);

	for (auto cur_pert : perturbation) {
		pert_A = A * (1 + cur_pert);
		pert_B = B * (1 + cur_pert);
		solution = SolveBoundaryProblem(a, b, n, alpha_0, alpha_1, pert_A, pert_B);
		outfile << cur_pert << " " << GetFactError(a, b, n, solution, f) << endl;
	}
	outfile.close();
}

int main(void) {
	double a = 0.1;
	//double b = 1;
	double b = 0.3;

	double alpha_0 = 1;
	double alpha_1 = 2;
	double A = alpha_0 * Solution(a) + alpha_1 * DerivativeSolution(a);
	double B = Solution(b);

	int n;

	n = 2;
	vector<double> res = SolveBoundaryProblem(a, b, n, alpha_0, alpha_1, A, B);
	vector<double> exact = SolutionToGridFunction(Solution, a, b, n);
	WriteGridFunctionToFile(a, b, n, "calc1.txt", res);
	WriteGridFunctionToFile(a, b, n, "sol1.txt", exact);

	n = 1024;
	res = SolveBoundaryProblem(a, b, n, alpha_0, alpha_1, A, B);
	exact = SolutionToGridFunction(Solution, a, b, n);
	WriteGridFunctionToFile(a, b, n, "calc2.txt", res);
	WriteGridFunctionToFile(a, b, n, "sol2.txt", exact);

	FormDependenceStepError(a, b, Solution, alpha_0, alpha_1, A, B);

	FormDependencePerturbationError(a, b, 10, Solution, alpha_0, alpha_1, A, B);

	return 0;
}
