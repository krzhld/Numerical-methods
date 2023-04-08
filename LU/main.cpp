#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

#define SIZE 3

typedef vector<double> column_t;
typedef vector<column_t> matrix_t;
typedef tuple<tuple<double, double, double>, matrix_t, column_t> system_t;
typedef pair<matrix_t, matrix_t> factorization_t;

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
			for (int k = 0; k < m ; k++)
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

int main(void) {
	matrix_t A(SIZE);
	A[0].push_back(6); A[0].push_back(3); A[0].push_back(45);
	A[1].push_back(12); A[1].push_back(13); A[1].push_back(94);
	A[2].push_back(18); A[2].push_back(37); A[2].push_back(159);

	//factorization_t fact = getLUFactorization(A);

	/*matrix_t res(SIZE);
	double temp;
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			temp = 0;
			for (int k = 0; k < SIZE; k++)
				temp += fact.first[i][k] * fact.second[k][j];
			res[i].push_back(temp);
		}
	}*/

	column_t ground_truth = { 1, 1, 1 };
	column_t column_B = { 54, 119, 214 };

	column_t sol = solveSystem(A, column_B);

	return 0;
}