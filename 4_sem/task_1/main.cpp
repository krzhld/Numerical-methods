#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double Function(double x) {
	return x * x * cos(2 * x) + 1;
}

pair<vector<double>, vector<double>> FormUniformGrid(int n, double a, double b) {
	vector<double> x(n + 1);
	vector<double> y(n + 1);
	double step = (b - a) / n;
	for (int i = 0; i <= n; i++) {
		x[i] = a + i * step;
		y[i] = Function(x[i]);
	}
	return make_pair(x, y);
}

vector<double> FormSeparatedDifferences(pair<vector<double>, vector<double>> grid) {
	vector<double> x = grid.first;
	vector<double> y = grid.second;
	int n = x.size();
	--n;
	vector<double> separatedDifferences(n + 1);
	vector<double> result(n + 1);
	for (int i = 0; i <= n; i++)
		separatedDifferences[i] = y[i];
	result[0] = y[0];
	for (int step = 1; step <= n; step++) {
		for (int i = 0; i <= n - step; i++)
			separatedDifferences[i] = (separatedDifferences[i + 1] - separatedDifferences[i]) / (x[i + step] - x[i]);
		result[step] = separatedDifferences[0];
	}
	separatedDifferences.clear();
	x.clear();
	y.clear();
	return result;
}

double NewtonPolynominal(int n, vector<double> separatedDifferences, vector<double> xNodes, double x) {
	double result = 0;
	double bracketComposition = 1;
	result += separatedDifferences[0];
	for (int step = 1; step <= n; step++) {
		bracketComposition *= (x - xNodes[step - 1]);
		result += separatedDifferences[step] * bracketComposition;
	}
	return result;
}

double FactError(double y1, double y2) {
	return fabs(y1 - y2);
}

vector<double> FactErrorInMiddlePoints(vector<double> separatedDifferences, vector<double> xNodes) {
	int n = xNodes.size();
	--n;
	double tempPoint;
	vector<double> resultErrors(n);
	for (int i = 0; i < n; i++) {
		tempPoint = (xNodes[i + 1] + xNodes[i]) / 2;
		resultErrors[i] = FactError(Function(tempPoint), NewtonPolynominal(n, separatedDifferences, xNodes, tempPoint));
	}
	return resultErrors;
}

bool FormInterpolatedPoints(int n, double step, double a, double b, string fileName) {
	ofstream outfile(fileName);
	if (!outfile)
		return false;
	pair<vector<double>, vector<double>> grid = FormUniformGrid(n, a, b);
	vector<double> separatedDifferences = FormSeparatedDifferences(grid);
	double temp, factError;
	for (double curPoint = a; curPoint <= b; curPoint += step) {
		temp = NewtonPolynominal(n, separatedDifferences, grid.first, curPoint);
		outfile << curPoint << " " << temp << endl;
	}
	grid.first.clear();
	grid.second.clear();
	separatedDifferences.clear();
	return true;
}

bool FormDependence(double step, double a, double b, string fileName) {
	ofstream outfile(fileName);
	if (!outfile)
		return false;
	pair<vector<double>, vector<double>> grid;
	vector<double> separatedDifferences;
	double maxError = 0, tempError;
	for (int n = 3; n <= 100; n++) {
		grid = FormUniformGrid(n, a, b);
		separatedDifferences = FormSeparatedDifferences(grid);
		for (double curPoint = a; curPoint <= b; curPoint += step) {
			tempError = FactError(NewtonPolynominal(n, separatedDifferences, grid.first, curPoint), Function(curPoint));
			if (tempError > maxError)
				maxError = tempError;
		}
		outfile << n << " " << maxError << endl;
		maxError = 0;
	}
	grid.first.clear();
	grid.second.clear();
	separatedDifferences.clear();
	return true;
}

int main(void) {
	/*FormInterpolatedPoints(3, 0.01, 0, 3, "data1.txt");
	FormInterpolatedPoints(5, 0.01, 0, 3, "data2.txt");
	FormInterpolatedPoints(7, 0.01, 0, 3, "data3.txt");*/
	FormDependence(0.05, 0, 3, "dependence.txt");
	return 0;
}
