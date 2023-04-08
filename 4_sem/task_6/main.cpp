#include <iostream>
#include <fstream>
#include <vector>

#define ORDER 3

using namespace std;

double function(double x, double y) {
	return (pow(x, 2) + pow(y, 3)) / (x * pow(y, 2));
	//return (y * y * exp(x) - 2 * y);
}

double solution(double x) {
	return pow((x - 1) * 3 * pow(x, 2), 1.0 / 3);
	//return exp(-x);
}

double getFactAccuracy(double(*f)(double), double a, double b, int n, vector<double>& grid_function) {
	double h = (b - a) / n;
	double result_accuracy = 0;
	double cur_accuracy = 0;

	for (int i = 0; i <= n; i++) {
		cur_accuracy = abs(grid_function[i] - f(a + i * h));
		if (cur_accuracy > result_accuracy)
			result_accuracy = cur_accuracy;
	}

	return result_accuracy;

}

// n - number of intervals
vector<double> predCorrScheme(double a, double b, double y_a, int n) {
	double h = (b - a) / n;
	vector<double> grid_function;
	vector<double> x;
	vector<double> y;
	vector<double> f;

	// getting start points with runge-kutta method 
	int i = 0;
	double x_i = a;
	x.push_back(x_i);
	double y_i = solution(a);
	y.push_back(y_i);
	grid_function.push_back(y_i);
	f.push_back(function(x_i, y_i));

	i++;
	for (; i < 3; i++) {
		double k1 = function(x_i, y_i);
		double k2 = function(x_i + h / 2, y_i + h * k1 / 2);
		double k3 = function(x_i + h, y_i - h * k1 + 2 * h * k2);

		y_i = y_i + h * (k1 + 4 * k2 + k3) / 6;
		grid_function.push_back(y_i);
		y.push_back(y_i);

		x_i = a + (double)i * h;
		x.push_back(x_i);

		f.push_back(function(x_i, y_i));
	}

	// pred-corr scheme
	double y_explicit;
	double y_implicit;
	for (; i <= n; i++) {
		y_explicit = (5 * f[0] - 16 * f[1] + 23 * f[2]) * h / 12 + y[2];
		y_implicit = (-f[1] + 8 * f[2] + 5 * function(x[2] + h, y_explicit)) * h / 12 + y[2];
		grid_function.push_back(y_implicit);

		y[0] = y[1]; 
		y[1] = y[2];
		y[2] = y_implicit;

		x[0] = x[1];
		x[1] = x[2];
		x[2] = a + i * h;

		f[0] = f[1];
		f[1] = f[2];
		f[2] = function(x[2], y[2]);
	}

	return grid_function;
}

vector<double> solutionToGridFunction(double(*f)(double), double a, double b, int n) {
	vector<double> grid_function;
	double h = (b - a) / n;

	for (int i = 0; i <= n; i++)
		grid_function.push_back(f(a + i * h));

	return grid_function;
}

void writeGridFunctionToFile(double a, double b, int n, string file_name, vector<double>& grid_function) {
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

void writeDependenceToFile(double(*f)(double), double a, double b, double y_a, string file_name) {
	ofstream outfile(file_name);
	outfile.precision(15);

	int number_steps = 18;
	int n = 3;
	double h;
	vector<double> grid_function;
	double cur_accuracy;
	double constant;
	while (number_steps--) {
		h = (b - a) / n;
		grid_function = predCorrScheme(a, b, y_a, n);
		cur_accuracy = getFactAccuracy(f, a, b, n, grid_function);
		constant = cur_accuracy / pow(h, ORDER);
		outfile << h << " " << cur_accuracy << " " << constant << endl;
		n *= 2;
	}
}

int main(void) {
	double a = 1.1; 
	double b = 3;
	double y_a = solution(a);

	vector<double> grid_function = predCorrScheme(a, b, y_a, 512);
	writeGridFunctionToFile(a, b, 512, "calc_sol_1.txt", grid_function);
	grid_function = solutionToGridFunction(solution, a, b, 512);
	writeGridFunctionToFile(a, b, 512, "exact_sol_1.txt", grid_function);

	grid_function = predCorrScheme(a, b, y_a, 1024);
	writeGridFunctionToFile(a, b, 1024, "calc_sol_2.txt", grid_function);
	grid_function = solutionToGridFunction(solution, a, b, 1024);
	writeGridFunctionToFile(a, b, 1024, "exact_sol_2.txt", grid_function);

	writeDependenceToFile(solution, a, b, y_a, "dependence.txt");

	return 0;
}
