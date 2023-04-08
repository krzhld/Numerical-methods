#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <vector>
#include <bits.h>

using namespace std;

inline double Function(double x) {
	return pow(x, 5) - 2.9 * pow(x, 3) + 6.5 * pow(x, 2) - 7 * x - 5.4;
}

double integral_truth = 30.789586666666667;

tuple<double, double, int, double> Integrate3Of8(double (*func)(double), double a, double b, double eps) {
	int n = 1;
	int powerN = 0;
	double h = (b - a) / (n * 3.0);

	double curIntegral = 3 * h / 8 * (func(a) + 3 * (func(a + h) + func(a + 2 * h)) + func(b));
	int p = 4;
	double denominatorRunge = (pow(2, p) - 1);
	double curPoint;
	double prevIntegral;
	do {
		prevIntegral = curIntegral;
		curIntegral = 0.0;
		n *= 2;
		++powerN;
		h = (b - a) / (n * 3.0);
		curPoint = a;
		for (int i = 0; i < n; i++) {
			curIntegral += (func(curPoint) + 3 * (func(curPoint + h) + func(curPoint + 2 * h)) + func(curPoint + 3 * h));
			curPoint += 3 * h;
		}
		curIntegral *= 3 * h / 8;
	} while (abs(curIntegral - prevIntegral) / denominatorRunge >= eps);

	return { curIntegral, abs(curIntegral - integral_truth), powerN, h };
}

double FindConstant(double (*func)(double), double a, double b) {
	// fourth derivative is 120 * x
	double derivative_constant = 120 * abs((a + b) / 2);
	double h = (b - a) / 3;
	double integral = (func(a) + 3 * (func(a + h) + func(a + 2 * h)) + func(b));
	integral *= 3 * h / 8;
	double fact_error = abs(integral_truth - integral);
	double constant = fact_error / pow(h, 5) / derivative_constant;
	return constant;
}

int main(void) {
	double a = -2.6; double b = -0.4;

	double constant;

	auto filename = "error_and_iter_on_eps.txt";
	ofstream outfile(filename);
	outfile.precision(15);
	auto eps = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
	double curNumIter, curError, curIntegral, curRank;
	for (auto curEps : eps) {
		tie(curIntegral, curError, curNumIter, curRank) = Integrate3Of8(Function, a, b, curEps);
		outfile << curEps << " " << curError << " " << curNumIter << " " << curRank << endl;
	}
	outfile.close();
	cout << "Constant: " << FindConstant(Function, a, b) << endl;

	return 0;
}
