#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

#define SIZE 10
#define NUMBER_EPSES 15

typedef vector<double> column_t;
typedef vector<column_t> matrix_t;
typedef tuple<tuple<double, double, double>, matrix_t, column_t> system_t;
typedef pair<matrix_t, matrix_t> factorization_t;

system_t readSystem(void) {
	FILE* fileID = fopen("MATLAB/matrix.txt", "rt");
	if (fileID == NULL) {
		cout << "Error reading file!" << endl;
		exit(-1);
	}

	/* lambda_1, lambda_2 and lambda_n */
	double lambda_1, lambda_2, lambda_n;
	fscanf(fileID, "%lf", &lambda_1);
	fscanf(fileID, "%lf", &lambda_2);
	fscanf(fileID, "%lf", &lambda_n);
	tuple<double, double, double> lambdas = make_tuple(lambda_1, lambda_2, lambda_n);

	double temp;
	matrix_t matrix_A;
	matrix_A.resize(SIZE);

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			fscanf(fileID, "%lf ", &temp);
			matrix_A[i].push_back(temp);
		}
	}

	column_t ground_truth;
	for (int i = 0; i < SIZE; i++) {
		fscanf(fileID, "%lf ", &temp);
		ground_truth.push_back(temp);
	}

	return make_tuple(lambdas, matrix_A, ground_truth);
}

double scalarComposition(column_t col1, column_t col2) {
	double res = 0;
	for (int i = 0; i < SIZE; i++)
		res += col1[i] * col2[i];
	return res;
}

double getSecNormColumn(column_t column) {
	return sqrt(scalarComposition(column, column));
}

matrix_t shiftMatrix(matrix_t matrix, double m) {
	for (int i = 0; i < SIZE; i++)
		matrix[i][i] -= m;
	return matrix;
}

column_t getOrthonormalVec(column_t col) {
	double norm = getSecNormColumn(col);
	column_t res(SIZE);
	if (norm == 0)
		return res;
	for (int i = 0; i < SIZE; i++)
		col[i] /= norm;
	return col;
}

pair<column_t, double> methodScalarComposition(matrix_t matrix, double eps, column_t first_approx, int* number_iter, double shift) {
	matrix_t B = shiftMatrix(matrix, shift);

	column_t x_k_1(SIZE);
	column_t x_k(SIZE);

	double lambda_k = 0, lambda_k_1;

	x_k = getOrthonormalVec(first_approx);

	*number_iter = 0;
	do {
		lambda_k_1 = lambda_k;
		x_k_1 = x_k;

		/* x_k = B * x_k_1 */
		for (int i = 0; i < SIZE; i++) {
			x_k[i] = 0;
			for (int j = 0; j < SIZE; j++)
				x_k[i] += B[i][j] * x_k_1[j];
		}

		lambda_k = scalarComposition(x_k, x_k) / scalarComposition(x_k, x_k_1);
		x_k = getOrthonormalVec(x_k);
		++(*number_iter);

	} while (abs(lambda_k - lambda_k_1) > eps);

	return { x_k, lambda_k + shift};
}

/*void printVector(column_t vector) {
	for (int i = 0; i < SIZE; i++)
		cout << vector[i] << " ";
	cout << endl;
}*/

double getSecNormOfDiscrepancy(matrix_t matrix, double eigen_value, column_t eigen_vector) {
	double norm = 0, temp = 0;
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++)
			temp += matrix[i][j] * eigen_vector[j];
		norm += pow(temp - eigen_value * eigen_vector[i], 2);
		temp = 0;
	}
	norm = sqrt(norm);
	return norm;
}

double getSecNormOfFactualError(column_t col1, column_t col2) {
	double norm = 0;
	for (int i = 0; i < SIZE; i++)
		norm += pow(col1[i] - col2[i], 2);
	norm = sqrt(norm);
	return norm;
}

void researchOptimalShift(matrix_t matrix) {
	column_t first_approx;
	srand(1);
	for (int i = 0; i < SIZE; i++)
		first_approx.push_back((double)rand() / RAND_MAX);

	double eps = 1e-6;
	int number_iter;

	pair<column_t, double> eigen_pair;

	FILE* fileID = fopen("MATLAB/research_optimal_shift.txt", "wt");
	if (fileID == NULL)
		return;

	for (double shift = 0.0; shift <= 60.0; shift += 5.0)
		fprintf(fileID, "%lf ", shift);
	fprintf(fileID, "\n");

	for (double shift = 0.0; shift <= 60.0; shift += 5.0) {
		eigen_pair = methodScalarComposition(matrix, eps, first_approx, &number_iter, shift);
		fprintf(fileID, "%d ", number_iter);
	}
	fclose(fileID);
}

factorization_t getLUFactorization(matrix_t matrix_A) {

	matrix_t matrix_L;
	matrix_L.resize(SIZE);
	for (int k = 0; k < SIZE; k++) {
		matrix_L[k].resize(SIZE);
		matrix_L[k][k] = 1;
	}

	matrix_t matrix_U;
	matrix_U.resize(SIZE);
	for (int k = 0; k < SIZE; k++)
		matrix_U[k].resize(SIZE);

	double temp_sum;

	//zeroth step

	for (int j = 0; j < SIZE; j++)
		matrix_U[0][j] = matrix_A[0][j];
	for (int i = 1; i < SIZE; i++)
		matrix_L[i][0] = matrix_A[i][0] / matrix_U[0][0];

	// other steps
	for (int m = 1; m < SIZE; m++) {
		for (int j = m - 1; j < SIZE; j++) {
			temp_sum = 0.0;
			for (int k = 0; k < m; k++)
				temp_sum += matrix_L[m][k] * matrix_U[k][j];
			matrix_U[m][j] = matrix_A[m][j] - temp_sum;
		}

		for (int i = m; i < SIZE; i++) {
			temp_sum = 0.0;
			for (int k = 0; k < m; k++)
				temp_sum += matrix_L[i][k] * matrix_U[k][m];
			matrix_L[i][m] = (matrix_A[i][m] - temp_sum) / matrix_U[m][m];
		}

	}
	return { matrix_L, matrix_U };
}

column_t solveLUSystem(factorization_t factorization, column_t column_B) {
	column_t temp_column;
	temp_column.resize(SIZE);
	column_t solution;
	solution.resize(SIZE);

	matrix_t matrix_L = factorization.first;
	matrix_t matrix_U = factorization.second;

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

	solution[SIZE - 1] = temp_column[SIZE - 1] / matrix_U[SIZE - 1][SIZE - 1];
	// other steps
	for (int i = SIZE - 2; i >= 0; i--) {
		for (int j = SIZE - 1; j > i; j--)
			temp_sum += matrix_U[i][j] * solution[j];

		solution[i] = (temp_column[i] - temp_sum) / matrix_U[i][i];
		temp_sum = 0;
	}

	return solution;
}

column_t solveSystem(matrix_t matrix_A, column_t column_B) {
	factorization_t factorization = getLUFactorization(matrix_A);
	column_t column_calc_solution = solveLUSystem(factorization, column_B);

	return column_calc_solution;
}

pair<column_t, double> reverseIterationsWithShift(matrix_t A, double start_lambda, column_t start_vec, double eps) {
	column_t x_k = getOrthonormalVec(start_vec);
	column_t x_k_1;
	column_t delta(SIZE);

	double lambda_k = start_lambda, lambda_k_1, temp;
	do {
		x_k_1 = x_k;
		lambda_k_1 = lambda_k;
		x_k = solveSystem(shiftMatrix(A, lambda_k_1), x_k_1);

		if (x_k[0] * x_k_1[0] < 0)
			for (int i = 0; i < SIZE; i++)
				x_k[i] *= -1.0;

		temp = 0;
		for (int i = 0; i < SIZE; i++)
			if (x_k[i] != 0)
				temp += x_k_1[i] / x_k[i];

		lambda_k = lambda_k_1 + temp / SIZE;

		x_k = getOrthonormalVec(x_k);

		for (int i = 0; i < SIZE; i++)
			delta[i] = x_k[i] - x_k_1[i];

	} while (getSecNormColumn(delta) > eps);

	return { x_k, lambda_k };
}

int main(void) {

	/*
	spectrum: 8, ..., 80, 93, 154 => shift: 0, 5, ..., 60; opt_shift = 50.5
	separability = 93 / 154 ~= 0.60
	*/

	system_t system = readSystem();
	tuple<double, double, double> lambdas;
	double lambda_1, lambda_2, lambda_n;
	matrix_t matrix_A;
	column_t ground_truth;

	tie(lambdas, matrix_A, ground_truth) = system;
	tie(lambda_1, lambda_2, lambda_n) = lambdas;
	
	//researchOptimalShift(matrix_A);

	double eps = 1e-6;
	int number_iter;
	column_t first_approx;
	first_approx.push_back(1);
	for (int i = 1; i < SIZE; i++)
		first_approx.push_back(0);
	
	pair<column_t, double> eig_pair;
	/*eig_pair = methodScalarComposition(matrix_A, eps, first_approx, &number_iter, (lambda_2 + lambda_n) / 2);
	eig_pair = reverseIterationsWithShift(matrix_A, eig_pair.second, eig_pair.first, eps);*/
	
	column_t epses = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15 };
	double norm_discrepancy, norm_error_vec, norm_error_value;

	FILE* fileID = fopen("MATLAB/dependences.txt", "a");
	if (fileID == NULL)
		return -1;
	for (int iter = 0; iter < NUMBER_EPSES; iter++) {
		eig_pair = methodScalarComposition(matrix_A, epses[iter], first_approx, &number_iter, (lambda_2 + lambda_n) / 2);
		eig_pair = reverseIterationsWithShift(matrix_A, eig_pair.second, eig_pair.first, epses[iter]);
		norm_discrepancy = getSecNormOfDiscrepancy(matrix_A, eig_pair.second, eig_pair.first);
		norm_error_vec = getSecNormOfFactualError(ground_truth, eig_pair.first);
		norm_error_value = abs(lambda_1 - eig_pair.second);
		fprintf(fileID, "%le %le %le %d %le\n", epses[iter], norm_error_vec, norm_discrepancy, number_iter, norm_error_value);
	}
	fclose(fileID);

	return 0;
}
