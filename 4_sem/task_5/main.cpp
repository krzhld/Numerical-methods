#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <float.h>

using namespace std;

#define RUNGE_DENOMINATOR 7.0 // (pow(2, 3) - 1) - 3-rd order of method

typedef pair<vector<double>, vector<double>> table_t;

double Function(double x, double y) {
	return (pow(x, 2) + pow(y, 3)) / (x * pow(y, 2));
}

double Solution(double x) {
	return pow((x - 1) * 3 * pow(x, 2), 1.0 / 3);
}

table_t SolutionToTableFunction(double (*solution)(double), double a, double b, int n) {
	table_t table_function;
	double h = (b - a) / n;
	double cur_point = a;
	for (int i = 0; i <= n; i++) {
		table_function.first.push_back(cur_point);
		table_function.second.push_back(solution(cur_point));
		cur_point += h;
	}
	return table_function;
}

bool WriteTableFunctionToFile(table_t table_function, string file_name) {
	ofstream outfile(file_name);
	outfile.precision(15);
	if (!outfile)
		return false;
	if (table_function.first.size() != table_function.second.size())
		return false;
	for (size_t i = 0; i < table_function.first.size(); i++)
		outfile << table_function.first[i] << " " << table_function.second[i] << endl;
	outfile.close();
	return true;
}

/*
method has 3-rd order
coef 1/2
*/
table_t RungeKuttaMethodStep(double a, double b, int n, double y_a) {
	vector<double> grid;
	vector<double> grid_function;
	double h = (b - a) / n;
	double k1, k2, k3;
	double x0 = a;
	double y0 = y_a, yi;
	for (int i = 0; i <= n; i++) {
		grid.push_back(x0);
		grid_function.push_back(y0);

		k1 = Function(x0, y0);
		k2 = Function(x0 + h / 2, y0 + h * k1 / 2);
		k3 = Function(x0 + h, y0 - h * k1 + 2 * h * k2);

		yi = y0 + h * (k1 + 4 * k2 + k3) / 6;

		y0 = yi;
		x0 = a + (double)(i + 1) * h;
	}

	return { grid, grid_function };
}

double GetFactError(table_t table_func) {
	vector<double> grid_x = table_func.first;
	vector<double> grid_y = table_func.second;
	double res = 0;
	double cur_error = 0;
	int size = grid_x.size();
	for (int i = 0; i < size; i++) {
		cur_error = abs(grid_y[i] - Solution(grid_x[i]));
		if (cur_error > res)
			res = cur_error;
	}
	return res;
}

tuple<table_t, int> RungeKuttaMethodEps(double a, double b, double y_a, int n, double eps) {
	table_t res_table;
	res_table.first.push_back(a);
	res_table.second.push_back(y_a);
	table_t temp_table;
	double l = a;
	double r = b;
	double y0 = y_a;
	int N = 0;
	int res_N = 0 ;

	double res1, res2;
	for (int k = 0; k < n; ++k) {
		r = a + (b - a) * (k + 1) / n;
		temp_table = RungeKuttaMethodStep(l, r, (int)pow(2, N), y0);
		res2 = temp_table.second[(int)pow(2, N)];
		/*temp_table = RungeKuttaMethodStep(l, r, 2 * (int)pow(2, N), y0);
		res2 = temp_table.second[2 * (int)pow(2, N)];*/
		//res1 = 9999999999.0;
		do {
			++N;
			res1 = res2;
			temp_table = RungeKuttaMethodStep(l, r, (int)pow(2, N), y0);
			res2 = temp_table.second[(int)pow(2, N)];
		} while (abs(res2 - res1) / RUNGE_DENOMINATOR >= eps);
		l = r;
		//y0 = res2 + (res2 - res1) / RUNGE_DENOMINATOR;
		y0 = res1;
		res_table.first.push_back(r);
		res_table.second.push_back(res1);
		//res_table.second.push_back(res2 + (res2 - res1) / RUNGE_DENOMINATOR);
		if (N > res_N)
			res_N = N;
		N = 0;
	}
	return { res_table, res_N };
}

int main(void) {
	double a = 1.1, b = 3;
	double y_a = Solution(a);

	// solve equations with two steps

	
	table_t table_function;

	table_function = RungeKuttaMethodStep(a, b, 16, y_a);
	WriteTableFunctionToFile(table_function, "step1_calculated_solution.txt");

	table_function = SolutionToTableFunction(Solution, a, b, 16);
	WriteTableFunctionToFile(table_function, "step1_exact_solution.txt");

	table_function = RungeKuttaMethodStep(a, b, 32, y_a);
	WriteTableFunctionToFile(table_function, "step2_calculated_solution.txt");

	table_function = SolutionToTableFunction(Solution, a, b, 32);
	WriteTableFunctionToFile(table_function, "step2_exact_solution.txt");
	

	// form dependences between eps and fact error, number iter

	
	ofstream outfile1("dependences.txt");
	table_t table;
	double fact_error;
	int number_iter;
	auto eps = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14 };
	for (auto cur_eps : eps) {
		tie(table, number_iter) = RungeKuttaMethodEps(a, b, y_a, 5, cur_eps);
		fact_error = GetFactError(table);
		outfile1 << cur_eps << " " << fact_error << " " << " " << number_iter << endl;
		table.first.clear();
		table.second.clear();
	}
	outfile1.close();
	

	// dependence between perturbation and fact error

	
	ofstream outfile2("perturbation.txt");
	table_t solution;
	int _;
	double pert_y_a;
	double fixed_eps = 1e-4;
	auto pert = { 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
	for (auto cur_pert : pert) {
		pert_y_a = y_a * (1.0 + cur_pert);
		tie(solution, _) = RungeKuttaMethodEps(a, b, pert_y_a, 5, fixed_eps);
		fact_error = GetFactError(solution);
		outfile2 << cur_pert << " " << fact_error << endl;
	}
	outfile2.close();
	

	//cout << Solution(1.3);
	return 0;
}
