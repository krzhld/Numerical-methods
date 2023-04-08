#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include <time.h>
using namespace std;

#define SIZE 15

typedef vector<double> column_t;
typedef vector<double> diag_matrix_t;
typedef vector<column_t> matrix_t;
typedef pair<matrix_t, diag_matrix_t> factorization_t;
typedef tuple<matrix_t, column_t, column_t, double> system_t;

// useless function
column_t multiplyMatrixColumn(matrix_t matrix, column_t column) {
	column_t answer;
	double temp;
	for (int i = 0; i < SIZE; i++) {
		temp = 0;
		for (int k = 0; k < SIZE; k++)
			temp += matrix[i][k] * column[k];
		answer.push_back(temp);
	}

	return answer;
}

system_t readFromFile(void) {
	FILE* fileID = fopen("MATLAB/matrixes.txt", "rt");
	if (fileID == NULL) {
		cout << "Error reading file!" << endl;
		exit(-1);
	}

	double cond_number;
	fscanf(fileID, "%lf", &cond_number);

	double temp;

	matrix_t matrix_A;
	matrix_A.resize(SIZE);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			fscanf(fileID, "%lf ", &temp);
			matrix_A[i].push_back(temp);
		}
	}

	column_t column_B;
	for (int i = 0; i < SIZE; i++) {
		fscanf(fileID, "%lf ", &temp);
		column_B.push_back(temp);
	}

	column_t column_ground_truth;
	for (int i = 0; i < SIZE; i++) {
		fscanf(fileID, "%lf ", &temp);
		column_ground_truth.push_back(temp);
	}

	fclose(fileID);

	system_t system = make_tuple(matrix_A, column_B, column_ground_truth, cond_number);
	return system;
}

void printMatrix(matrix_t matrix) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void printVector(vector<double> vec) {
	for (int j = 0; j < SIZE; j++) {
		cout << vec[j] << " ";
	}
	cout << endl;
}

factorization_t getLDLtFactorization(matrix_t matrix_A) {

	matrix_t matrix_L;
	matrix_L.resize(SIZE);
	for (int k = 0; k < SIZE; k++) {
		matrix_L[k].resize(SIZE);
		matrix_L[k][k] = 1;
	}
		
	diag_matrix_t matrix_D;
	matrix_D.resize(SIZE);

	// zeroth step
	matrix_D[0] = matrix_A[0][0];
	for (int i = 1; i < SIZE; i++)
		matrix_L[i][0] = matrix_A[i][0] / matrix_D[0];

	double temp_sum;

	// other steps
	for (int m = 1; m < SIZE; m++) {
		temp_sum = 0;
		// calculate sum for diag element
		for (int k = 0; k < m; k++)
			temp_sum += matrix_L[m][k] * matrix_L[m][k] * matrix_D[k];

		matrix_D[m] = matrix_A[m][m] - temp_sum;
		
		// calculate m-column of L matrix
		for (int i = m; i < SIZE; i++) {
			temp_sum = 0;
			// calculate auxiliary sum 
			for (int k = 0; k < m; k++) 
				temp_sum += matrix_L[i][k] * matrix_D[k] * matrix_L[m][k];
			if (i != m)
				matrix_L[i][m] = (matrix_A[i][m] - temp_sum) / matrix_D[m];
		}

	}
	return { matrix_L, matrix_D };
}

column_t solveLDLtSystem(factorization_t factorization, column_t column_B) {
	column_t temp_column;
	temp_column.resize(SIZE);
	column_t solution;
	solution.resize(SIZE);

	matrix_t matrix_L = factorization.first;
	column_t matrix_D = factorization.second;

	double temp_sum = 0;

	/* 1 epoch */
	temp_column[0] = column_B[0];

	// other steps
	for (int i = 1; i < SIZE; i++) {
		// calculate auxiliary sum
		for (int j = 0; j < i; j++)
			temp_sum += matrix_L[i][j] * temp_column[j];

		temp_column[i] = (column_B[i] - temp_sum);
		temp_sum = 0;
	}

	/* 2 epoch */
	for (int i = 0; i < SIZE; i++)
		temp_column[i] = temp_column[i] / matrix_D[i];

	/* 3 epoch */
	solution[SIZE - 1] = temp_column[SIZE - 1];

	// other steps
	for (int i = SIZE - 2; i >= 0; i--) {
		// calculate auxiliary sum
		for (int j = SIZE - 1; j > i; j--)
			temp_sum += matrix_L[j][i] * solution[j];

		solution[i] = temp_column[i] - temp_sum;
		temp_sum = 0;
	}

	return solution;
}

column_t solveSystem(system_t system) {
	matrix_t matrix_A; column_t column_B, _; double cond_number;
	tie(matrix_A, column_B, _, cond_number) = system;

	factorization_t factorization = getLDLtFactorization(matrix_A);
	column_t column_calc_solution = solveLDLtSystem(factorization, column_B);

	return column_calc_solution;
}

double norm2(column_t column) {
	double norm = 0;
	for (int i = 0; i < SIZE; i++) {
		norm += pow(column[i], 2);
	}
	norm = sqrt(norm);
	return norm;
}

double getFactError(column_t ground_truth_solution, column_t calc_solution, bool print_error) {
	double error = 0;
	for (int i = 0; i < SIZE; i++) {
		error += pow(ground_truth_solution[i] - calc_solution[i], 2);
	}
	error = sqrt(error);

	if (print_error)
		cout << "Actual error: " << error << endl;

	return error;
}

double getDiscrepancy(matrix_t matrix_A, column_t calc_solution, column_t column_B, bool print_discrepancy) {
	double discrepancy = 0;
	column_t left_part_of_equation = multiplyMatrixColumn(matrix_A, calc_solution);
	for (int i = 0; i < SIZE; i++) {
		discrepancy += pow(left_part_of_equation[i] - column_B[i], 2);
	}
	discrepancy = sqrt(discrepancy);

	if (print_discrepancy)
		cout << "Discrepancy: " << discrepancy << endl;

	return discrepancy;
}

column_t getPerturbedColumn(column_t column) {
	column_t perturbed_column(SIZE);
	double coefficient = 0;
	srand(time(NULL));
	for (int i = 0; i < SIZE; i++) {
		coefficient = ((rand() % 11) - 5) / (100 * 1.0);
		perturbed_column[i] = (column[i] * (1 + coefficient));
	}
	return perturbed_column;
}

bool checkInequality(system_t system, column_t perturbedRightPart, column_t calc_solution) {
	matrix_t matrix_A; column_t column_B, ground_truth, _; double cond_number = 0e1;
	tie(matrix_A, column_B, ground_truth, cond_number) = system;

	system_t perturbed_system = make_tuple(matrix_A, perturbedRightPart, _, cond_number);
	column_t perturber_calc_solution = solveSystem(perturbed_system);

	double dx = 0.0, x = 0.0, db = 0.0, b = 0.0;

	for (int i = 0; i < SIZE; i++) {
		dx += pow(perturber_calc_solution[i] - calc_solution[i], 2);
		x += pow(calc_solution[i], 2);
			db += pow(perturbedRightPart[i] - column_B[i], 2);
		b += pow(column_B[i], 2);
	}
	dx = sqrt(dx);
	x = sqrt(x);
	db = sqrt(db);
	b = sqrt(b);

	if (dx / x <= cond_number * db / b) {
		cout << "dx: " << dx << endl << "x: " << x << endl << "db: " << db << endl << "b: " << b << endl;
		return true;
	}
	else
		return false;
}

void writeDependence(double cond_number, double fact_error, double discrepancy) {
	FILE* fileID = fopen("MATLAB/dependence.txt", "a");
	if (fileID == NULL) {
		cout << "Error reading file!" << endl;
		exit(-1);
	}

	fprintf(fileID, "%le %le %le\n", cond_number, fact_error, discrepancy);
}

int main() {
	/* system with certain condition number */
	system_t system = readFromFile();
	column_t calc_solution = solveSystem(system);

	matrix_t matrix_A; column_t column_B, ground_truth; double cond_number;
	tie(matrix_A, column_B, ground_truth, cond_number) = system;

	printVector(ground_truth);
	printVector(calc_solution);

	double fact_error = getFactError(ground_truth, calc_solution, true);
	double discrepancy = getDiscrepancy(matrix_A, calc_solution, column_B, true);

	if (checkInequality(system, getPerturbedColumn(column_B), calc_solution))
		cout << "true" << endl;
	else
		cout << "false" << endl;

	//writeDependence(cond_number, fact_error, discrepancy);

	return 0;
}
