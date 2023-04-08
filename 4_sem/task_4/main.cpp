#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

using namespace std;

double groundTruth = 23.83741577002152;

double Function(double x) {
	return cos(0.4 * x) * (pow(x, 5) - 2.9 * pow(x, 3) + 6.5 * pow(x, 2) - 7 * x - 5.4);
}

/*double Primitive(double x) {
	double res;
	res = (31.25 * pow(x, 4) - 2398.12 * pow(x, 2) + 81.25 * x + 29932.8) * cos(0.4 * x);
	res += (2.5 * pow(x, 5) - 319.75 * pow(x, 3) + 16.25 * pow(x, 2) + 11973.1 * x - 216.625) * sin(0.4 * x);
	return res;
}*/

tuple<double, double, double> CalculateCoeff(double a, double b) {
	double B1 = b - a;
	double discriminant = abs(b - a) * sqrt(3) / 6;
	double c = (a + b) / 2 + discriminant;
	double B = c - (b + a) / 2;
	return { B, c, B1 };
}

tuple<double, int, double> IntegrateRalston(double (*func)(double), double a, double b, double eps) {
	int powerN = 0;
	int N = 0;
	double h = (b - a);
	double B, x1, B1;
	tie(B, x1, B1) = CalculateCoeff(a, b);
	double curIntegral = B * (Function(a) - Function(b)) + B1 * Function(x1);
	double prevIntegral = 0;
	double curPoint = 0;
	double denominatorRunge = (pow(2, 3) - 1); // 3-rd order of method
	do {
		prevIntegral = curIntegral;
		curIntegral = 0.0;
		++powerN;
		N = pow(2, powerN);
		h = (b - a) / N;
		curPoint = a;
		for (int i = 0; i < N; i++) {
			tie(B, x1, B1) = CalculateCoeff(curPoint, curPoint + h);
			curIntegral += B * (Function(curPoint) - Function(curPoint + h)) + B1 * Function(x1);
			curPoint += h;
		}
	} while (abs(curIntegral - prevIntegral) / denominatorRunge >= eps);

	// Richardson Correction
	double correctedIntegral = curIntegral + (curIntegral - prevIntegral) / denominatorRunge;

	return { curIntegral, powerN, correctedIntegral };
}


int main(void) {
	double a = -2.6;
	double b = -0.4;
	auto filename = "error_and_iter_on_eps.txt";
	ofstream outfile(filename);
	outfile.precision(12);
	auto eps = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
	double curNumIter, curIntegral, correctedIntegral;
	for (auto curEps : eps) {
		tie(curIntegral, curNumIter, correctedIntegral) = IntegrateRalston(Function, a, b, curEps);
		outfile << curEps << " " << abs(curIntegral - groundTruth) << " " << curNumIter << " " << abs(correctedIntegral - groundTruth) << endl;
	}
	outfile.close();
	return 0;
}
