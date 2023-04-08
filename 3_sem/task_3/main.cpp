#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <iostream>
#include <tuple>
#include <cmath>

using namespace std;

#define SIZE 15
#define NUMBER_EPSES 15

typedef vector<double> column_t;
typedef vector<column_t> matrix_t;
typedef tuple<matrix_t, column_t, column_t> SOLE_t;
typedef pair<vector<SOLE_t>, double> systems_t;

systems_t readSystems(void) {
	FILE* fileID = fopen("MATLAB/SOLE.txt", "rt");
	if (fileID == NULL) {
		cout << "Error reading file!" << endl;
		exit(-1);
	}

	double det;
	fscanf(fileID, "%lf", &det);

	double temp;
	matrix_t matrix_A;
	matrix_A.resize(SIZE);
	column_t column_ground_truth, column_B;

	vector<SOLE_t> systems;

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			fscanf(fileID, "%lf ", &temp);
			matrix_A[i].push_back(temp);
		}
	}

	for (int i = 0; i < SIZE; i++) {
		fscanf(fileID, "%lf ", &temp);
		column_ground_truth.push_back(temp);
	}

	for (int i = 0; i < SIZE; i++) {
		fscanf(fileID, "%lf ", &temp);
		column_B.push_back(temp);
	}

	systems.push_back(make_tuple(matrix_A, column_ground_truth, column_B));

	for (int iter = 1; iter < NUMBER_EPSES; iter++) {
		for (int i = 0; i < SIZE; i++) {
			for (int j = 0; j < SIZE; j++) {
				fscanf(fileID, "%lf ", &temp);
				matrix_A[i][j] = temp;
			}
		}

		for (int i = 0; i < SIZE; i++) {
			fscanf(fileID, "%lf ", &temp);
			column_ground_truth[i] = temp;
		}

		for (int i = 0; i < SIZE; i++) {
			fscanf(fileID, "%lf ", &temp);
			column_B[i] = temp;
		}

		systems.push_back(make_tuple(matrix_A, column_ground_truth, column_B));
	}
	fclose(fileID);

	return {systems, det};
}

double getInfNormMatrix(matrix_t matrix) {
	double res = 0;
	double norm_row = 0;
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++)
			norm_row += abs(matrix[i][j]);
		if (norm_row > res)
			res = norm_row;
		norm_row = 0;
	}
	return res;
}

double getInfNormColumn(column_t column) {
	double res = 0;
	for (int i = 0; i < SIZE; i++)
		if (abs(column[i]) > res)
			res = abs(column[i]);
	return res;
}

pair<matrix_t, column_t> prepareSeidelMethod(matrix_t matrix_A, column_t column_B) {
	matrix_t C;
	C.resize(SIZE);
	double norm_A = getInfNormMatrix(matrix_A);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < i; j++)
			C[i].push_back(- matrix_A[i][j] / norm_A);

		C[i].push_back(1 - matrix_A[i][i] / norm_A);

		for (int j = i + 1; j < SIZE; j++)
			C[i].push_back(-matrix_A[i][j] / norm_A);
	}

	column_t g;
	for (int i = 0; i < SIZE; i++)
		g.push_back(column_B[i] / norm_A);
	
	return { C, g };
}

void printColumn(column_t column) {
	for (int i = 0; i < SIZE; i++)
		cout << column[i] << " ";
	cout << endl;
}

double getFactError(column_t ground_truth_solution, column_t calc_solution, bool print_error) {
	column_t fact_error (SIZE);
	for (int i = 0; i < SIZE; i++)
		fact_error[i] = calc_solution[i] - ground_truth_solution[i];

	double norm_fact_error = getInfNormColumn(fact_error);
	if (print_error)
		cout << "Actual error: " << norm_fact_error << endl;

	return norm_fact_error;
}

double getDiscrepancy(matrix_t matrix_A, column_t calc_solution, column_t column_B, bool print_discrepancy) {
	column_t left_part_of_equation (SIZE);

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
			left_part_of_equation[i] += matrix_A[i][j] * calc_solution[j];

	column_t discrepancy (SIZE);
	for (int i = 0; i < SIZE; i++)
		discrepancy[i] = left_part_of_equation[i] - column_B[i];

	double norm_discrepancy = getInfNormColumn(discrepancy);

	if (print_discrepancy)
		cout << "Discrepancy: " << norm_discrepancy << endl;

	return norm_discrepancy;
}

void writeDependence(double eps, double fact_error, double discrepancy, int number_iter) {
	FILE* fileID = fopen("MATLAB/dependence_edited.txt", "a");
	if (fileID == NULL) {
		cout << "Error reading file!" << endl;
		exit(-1);
	}

	fprintf(fileID, "%le %le %le %d\n", eps, fact_error, discrepancy, number_iter);
}

/* inf-norm */	
bool converge(column_t x_k, column_t x_k_1, double norm_C, double eps) {
	column_t delta;
	for (int i = 0; i < SIZE; i++)
		delta.push_back(x_k_1[i] - x_k[i]);
	double norm_delta = getInfNormColumn(delta);
	if (norm_C <= 0.5)
		return (norm_delta < eps);
	else
		return (norm_delta < (1 - norm_C) / norm_C * eps);
}

column_t seidelMethod(matrix_t C, column_t g, column_t firstApprox, double eps, int* number_iteration, bool printCalculations) {
	column_t x = firstApprox;
	column_t x_before;

	double norm_C = getInfNormMatrix(C);
	double temp = 0;
	int epoches = 0;
	do {
		x_before = x;
		/*if (printCalculations)
			printColumn(x);*/
		for (int i = 0; i < SIZE; i++) {
			for (int j = 0; j < SIZE; j++)
				temp += C[i][j] * x[j];
			temp += g[i];
			x[i] = temp;
			temp = 0;
		}
		++epoches;
	} while (!converge(x_before, x, norm_C, eps));
	if (number_iteration != NULL)
		*number_iteration = epoches;
	return x;
}

int main(void) {
	int epoches;

	column_t initial_approximation (SIZE);
	column_t x, column_B, ground_truth, g;
	column_t epses = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15 };

	matrix_t matrix_A, C;

	systems_t systems = readSystems();
	vector<SOLE_t> SOLEs = systems.first;
	double det = systems.second;

	pair<matrix_t, column_t> utility_elements;

	for (int iter = 0; iter < NUMBER_EPSES; iter++) {
		tie(matrix_A, ground_truth, column_B) = SOLEs[iter];
		utility_elements = prepareSeidelMethod(matrix_A, column_B);
		C = utility_elements.first;
		g = utility_elements.second;
		x = initial_approximation;
		epoches = 0;
		x = seidelMethod(C, g, x, epses[iter], &epoches, false);
		writeDependence(epses[iter], getFactError(ground_truth, x, false), getDiscrepancy(matrix_A, x, column_B, false), epoches);
		printColumn(x);
	}
	
	return 0;
}
