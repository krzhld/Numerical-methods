#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double Function(double x) {
	return x * x * cos(2 * x) + 1;
}

/*double DerivativeOfFunction(double x) {
	return 2 * x * (cos(2 * x) - x * sin(2 * x));
}*/

double FactError(double y1, double y2) {
	return abs(y1 - y2);
}

pair<vector<double>, vector<double>> FormUniformTable(int n, double a, double b) {
	vector<double> x(n + 1);
	vector<double> y(n + 1);
	double step = (b - a) / n;
	for (int i = 0; i <= n; i++) {
		x[i] = a + i * step;
		y[i] = Function(x[i]);
	}
	return make_pair(x, y);
}

vector<double> GEM(vector<vector<double>> matrix, vector<double> column) {
	int n = column.size();
	double maxElem = 0; int indexMaxElem = 0;
	vector<double> answer(n);
	vector<double> tempVec(n); 
	double temp;
	for (int step = 0; step < n; step++) {
		// find max elem in column
		for (int i = step; i < n; i++) {
			if (abs(matrix[i][step]) > maxElem) {
				maxElem = matrix[i][step];
				indexMaxElem = i;
			}
		}
		// swap rows
		if (maxElem != 0 && indexMaxElem != step) {
			tempVec = matrix[step];
			matrix[step] = matrix[indexMaxElem];
			matrix[indexMaxElem] = tempVec;
			temp = column[step];
			column[step] = column[indexMaxElem];
			column[indexMaxElem] = temp;
		}
		maxElem = 0;
		indexMaxElem = 0;
		// getting zeros under diag elem
		for (int i = step + 1; i < n; i++) {
			temp = matrix[i][step] / matrix[step][step];
			for (int j = step; j < n; j++)
				matrix[i][j] = matrix[i][j] - temp * matrix[step][j];
			column[i] = column[i] - temp * column[step];
		}
	}
	// reverse gear
	for (int i = n - 1; i >= 0; i--) {
		temp = 0;
		for (int j = i + 1; j < n; j++)
			temp += matrix[i][j] * answer[j];
		answer[i] = (column[i] - temp) / matrix[i][i];
	}
	tempVec.clear();
	return answer;
}

// n >= 3
vector<vector<double>> FormNotAKnotQuadraticSpline(pair<vector<double>, vector<double>> table) {
	int n = table.first.size();
	vector<double> grid = table.first;
	vector<double> gridFunc = table.second;
	vector<vector<double>> res(n - 2);

	vector<double> bufferColumn(3);
	vector<vector<double>> bufferMatrix(3);
	vector<double> bufferAnswer(3);

	// iteration with lost knot
	bufferMatrix = { { pow(grid[0], 2), grid[0], 1 }, { pow(grid[1], 2), grid[1], 1 }, { pow(grid[2], 2), grid[2], 1 } };
	bufferColumn = { gridFunc[0], gridFunc[1], gridFunc[2] };
	bufferAnswer = GEM(bufferMatrix, bufferColumn);
	res[0] = { grid[0], grid[2], bufferAnswer[0], bufferAnswer[1], bufferAnswer[2] };

	// other iterations
	// maybe problem with index in grid
	for (int i = 1; i < n - 2; i++) {
		bufferMatrix = { { pow(grid[i + 1], 2), grid[i + 1], 1 }, { pow(grid[i + 2], 2), grid[i + 2], 1 }, { 2 * grid[i + 1], 1, 0 } };
		bufferColumn = { gridFunc[i + 1], gridFunc[i + 2], (2 * bufferAnswer[0] * grid[i + 1] + bufferAnswer[1]) };
		bufferAnswer = GEM(bufferMatrix, bufferColumn);
		res[i] = { grid[i + 1], grid[i + 2], bufferAnswer[0], bufferAnswer[1], bufferAnswer[2] };
	}

	grid.clear();
	gridFunc.clear();
	bufferAnswer.clear();
	bufferColumn.clear();
	bufferMatrix.clear();

	return res;
}

double SplineFunction(double x, vector<vector<double>> coeff) {
	size_t iter = 0;
	while (iter < coeff.size()) {
		if (coeff[iter][0] < x && x < coeff[iter][1])
			return coeff[iter][2] * pow(x, 2) + coeff[iter][3] * x + coeff[iter][4];
		else
			++iter;
	}
	return (double)NAN;
}

bool FormInterpolatedPoints(vector<vector<double>> coeff, double step, string fileName) {
	ofstream outfile(fileName);
	if (!outfile)
		return false;
	double temp = 0;
	int iter = 0;
	double a = coeff[0][0];
	double b = coeff[coeff.size() - 1][1];
	for (double curPoint = a; curPoint <= b; curPoint += step) {
		if (curPoint > coeff[iter][1])
			++iter;
		temp = coeff[iter][2] * pow(curPoint, 2) + coeff[iter][3] * curPoint + coeff[iter][4];
		outfile << curPoint << " " << temp << endl;
	}
	return true;
}

bool FormDependence(double step, double a, double b, string fileName) {
	ofstream outfile(fileName);
	if (!outfile)
		return false;
	pair<vector<double>, vector<double>> table;
	vector<vector<double>> coeff;
	double maxError = 0, tempError, temp;
	int iter;
	for (int n = 3; n <= 100; n++) {
		iter = 0;
		table = FormUniformTable(n, a, b);
		coeff = FormNotAKnotQuadraticSpline(table);
		for (double curPoint = a; curPoint <= b; curPoint += step) {
			if (curPoint > coeff[iter][1])
				++iter;
			temp = coeff[iter][2] * pow(curPoint, 2) + coeff[iter][3] * curPoint + coeff[iter][4];
			tempError = FactError(temp, Function(curPoint));
			if (tempError > maxError)
				maxError = tempError;
		}
		outfile << n << " " << maxError << endl;
		maxError = 0;
		iter = 0;
		table.first.clear(); table.second.clear();
		coeff.clear();
	}
	return true;
}

int main(void) {
	//pair<vector<double>, vector<double>> table;
	//vector<vector<double>> coeff;

	/*table = FormUniformTable(3, 0, 3);
	coeff = FormNotAKnotQuadraticSpline(table);
	FormInterpolatedPoints(coeff, 0.01, "data1.txt");
	table.first.clear(); table.second.clear();
	coeff.clear();

	table = FormUniformTable(5, 0, 3);
	coeff = FormNotAKnotQuadraticSpline(table);
	FormInterpolatedPoints(coeff, 0.01, "data2.txt");
	table.first.clear(); table.second.clear();
	coeff.clear();

	table = FormUniformTable(7, 0, 3);
	coeff = FormNotAKnotQuadraticSpline(table);
	FormInterpolatedPoints(coeff, 0.01, "data3.txt");
	table.first.clear(); table.second.clear();
	coeff.clear();*/

	//FormDependence(0.01, 0.0, 3.0, "dependence.txt");

	return 0;
}
