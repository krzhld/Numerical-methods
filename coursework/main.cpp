#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>

using namespace std;

double function(double x, double y) {
	//return (pow(x, 2) + pow(y, 3)) / (x * pow(y, 2));
	return (x / y / y + y / x);
	//return -(1 + x * y) / (x * x);
	//return (y * y * exp(x) - 2 * y);
}

double solution(double x) {
	return pow((x - 1) * 3 * pow(x, 2), 1.0 / 3);
	//return (1 - log(abs(x))) / x;
	//return exp(-x);
}

tuple<double, double> getInfNormInaccuracy(double(*f)(double), double a, double b, int n, vector<double>& grid_function) {
	double h = (b - a) / n;
	double result_inaccuracy = 0;
	double cur_inaccuracy = 0;
	double cur_coord = a;
	double max_coord = a;

	for (int i = 0; i <= n; i++) {
		cur_inaccuracy = abs(grid_function[i] - f(cur_coord));
		if (cur_inaccuracy > result_inaccuracy) {
			max_coord = cur_coord;
			result_inaccuracy = cur_inaccuracy;
		}
		cur_coord += h;
	}

	return {result_inaccuracy, max_coord};
}

vector<double> rungeKuttaMethodStep(double(*f)(double, double), double a, double b, int n, double y_a) {
	vector<double> grid_function;
	double h = (b - a) / n;
	double k1, k2, k3;
	double x0 = a;
	double y0 = y_a, yi;
	for (int i = 0; i <= n; i++) {
		grid_function.push_back(y0);

		k1 = f(x0, y0);
		k2 = f(x0 + h / 2, y0 + h * k1 / 2);
		k3 = f(x0 + h, y0 - h * k1 + 2 * h * k2);

		yi = y0 + h * (k1 + 4 * k2 + k3) / 6;

		y0 = yi;
		x0 += h;
	}

	return grid_function;
}

vector<double> predCorrScheme(double(*f)(double, double), double a, double b, int n, double y_a) {
	double h = (b - a) / n;
	vector<double> grid_function;
	vector<double> x;
	vector<double> y;
	vector<double> func;

	// getting start points with runge-kutta method 
	int i = 0;
	double x_i = a;
	x.push_back(x_i);
	double y_i = solution(a);
	y.push_back(y_i);
	grid_function.push_back(y_i);
	func.push_back(f(x_i, y_i));

	i++;
	for (; i < 3; i++) {
		double k1 = f(x_i, y_i);
		double k2 = f(x_i + h / 2, y_i + h * k1 / 2);
		double k3 = f(x_i + h, y_i - h * k1 + 2 * h * k2);

		y_i = y_i + h * (k1 + 4 * k2 + k3) / 6;
		grid_function.push_back(y_i);
		y.push_back(y_i);

		x_i = a + (double)i * h;
		x.push_back(x_i);

		func.push_back(f(x_i, y_i));
	}

	// pred-corr scheme
	double y_explicit;
	double y_implicit;
	for (; i <= n; i++) {
		y_explicit = (5 * func[0] - 16 * func[1] + 23 * func[2]) * h / 12 + y[2];
		y_implicit = (-func[1] + 8 * func[2] + 5 * f(x[2] + h, y_explicit)) * h / 12 + y[2];
		grid_function.push_back(y_implicit);

		y[0] = y[1];
		y[1] = y[2];
		y[2] = y_implicit;

		x[0] = x[1];
		x[1] = x[2];
		x[2] = a + i * h;

		func[0] = func[1];
		func[1] = func[2];
		func[2] = f(x[2], y[2]);
	}

	return grid_function;
	/*double h = (b - a) / n;
	vector<double> grid_function;
	vector<double> x;
	vector<double> y;

	// getting start points with runge-kutta method 
	int i = 0;
	double x_i = a;
	x.push_back(x_i);
	double y_i = solution(a);
	y.push_back(y_i);
	grid_function.push_back(y_i);

	i++;
	for (; i < 4; i++) {
		double k1 = f(x_i, y_i);
		double k2 = f(x_i + h / 2, y_i + h * k1 / 2);
		double k3 = f(x_i + h, y_i - h * k1 + 2 * h * k2);

		y_i = y_i + h * (k1 + 4 * k2 + k3) / 6;
		grid_function.push_back(y_i);
		y.push_back(y_i);

		x_i = a + (double)i * h;
		x.push_back(x_i);
	}

	// pred-corr scheme
	double y_explicit;
	double y_implicit;
	for (; i <= n; i++) {
		y_explicit = (-9 * f(x[0], y[0]) + 37 * f(x[1], y[1]) - 59 * f(x[2], y[2]) + 55 * f(x[3], y[3])) * h / 24 + y[3];
		y_implicit = (9 * f(x[3] + h, y_explicit) + 19 * f(x[3], y[3]) - 5 * f(x[2], y[2]) + f(x[1], y[1])) * h / 24 + y[3];
		grid_function.push_back(y_implicit);

		y[0] = y[1];
		y[1] = y[2];
		y[2] = y[3];
		y[3] = y_implicit;

		x[0] = x[1];
		x[1] = x[2];
		x[2] = x[3];
		x[3] = a + (double)i * h;
	}
	x.clear();
	y.clear();

	return grid_function;*/
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

void writeErrorsToFile(double a, double b, int n, string file_name, vector<double>& grid_function_sol, vector<double>& grid_function_calc) {
	double h = (b - a) / n;
	double x = a;

	ofstream outfile(file_name);
	outfile.precision(15);
	for (int i = 0; i <= n; i++) {
		x = a + i * h;
		outfile << x << " " << (grid_function_calc[i] - grid_function_sol[i]) << endl;
	}
	outfile.close();
}

int main(void) {
	double a = 1.1;
	double b = 3;
	double y_a = solution(a);
	vector<double> grid_function_sol;
	vector<double> grid_function_calc;
	int n;


	n = 16;
	grid_function_sol = solutionToGridFunction(solution, a, b, n);
	writeGridFunctionToFile(a, b, n, "sol_exact_1.txt", grid_function_sol);
	grid_function_calc = rungeKuttaMethodStep(function, a, b, n, y_a);
	writeGridFunctionToFile(a, b, n, "sol_rungekutta_1.txt", grid_function_calc);
	writeErrorsToFile(a, b, n, "errors_rungekutta_1.txt", grid_function_sol, grid_function_calc);
	grid_function_calc = predCorrScheme(function, a, b, n, y_a);
	writeGridFunctionToFile(a, b, n, "sol_predcorr_1.txt", grid_function_calc);
	writeErrorsToFile(a, b, n, "errors_predcorr_1.txt", grid_function_sol, grid_function_calc);
	
	
	n = 32;
	grid_function_sol = solutionToGridFunction(solution, a, b, n);
	writeGridFunctionToFile(a, b, n, "sol_exact_2.txt", grid_function_sol);
	grid_function_calc = rungeKuttaMethodStep(function, a, b, n, y_a);
	writeGridFunctionToFile(a, b, n, "sol_rungekutta_2.txt", grid_function_calc);
	writeErrorsToFile(a, b, n, "errors_rungekutta_2.txt", grid_function_sol, grid_function_calc);
	grid_function_calc = predCorrScheme(function, a, b, n, y_a);
	writeGridFunctionToFile(a, b, n, "sol_predcorr_2.txt", grid_function_calc);
	writeErrorsToFile(a, b, n, "errors_predcorr_2.txt", grid_function_sol, grid_function_calc);

	grid_function_sol.clear();

	ofstream outfile("dependences.txt");
	outfile.precision(15);
	double inaccuracy_predcorr, inaccuracy_rungekutta;
	double coord_predcorr, coord_rungekutta;
	for (int i = 4; i <= 131072; i *= 2) {
		grid_function_calc = rungeKuttaMethodStep(function, a, b, i, y_a);
		tie(inaccuracy_rungekutta, coord_rungekutta) = getInfNormInaccuracy(solution, a, b, i, grid_function_calc);
		grid_function_calc = predCorrScheme(function, a, b, i, y_a);
		tie(inaccuracy_predcorr, coord_predcorr) = getInfNormInaccuracy(solution, a, b, i, grid_function_calc);
		outfile << i << " " << inaccuracy_rungekutta << " " << inaccuracy_predcorr << " " << coord_rungekutta << " " << coord_predcorr << endl;
	}
	outfile.close();
	grid_function_calc.clear();

	return 0;
}
