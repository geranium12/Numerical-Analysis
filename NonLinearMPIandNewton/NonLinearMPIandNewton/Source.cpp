#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double eps7 = 1E-7;
const double eps2 = 1E-2;

double fx(double x) {
	return pow(2, x) - x * x - 0.5;
}

double dfx(double x) {
	return pow(2, x) * log(2) - 2 * x;
}

double gx(double x) {
	return -sqrt(pow(2, x) - 0.5);
}

typedef double(*function) (double x);

double dichotomy(function fx, double a, double b) {
	cout << "DICHOTOMY" << endl;
	int k = 0;
	double m = (a + b) / 2;
	cout << fixed << setprecision(6);
	cout << "k     " << "a             " << "b            " << "(a + b) / 2   " << "f(a)          " << "f(b)         " << "f((a + b) / 2)     " << "b - a     " << endl;
	while (fabs((b - a) / 2) >= eps2) {
		cout << k << "     " << a << "     " << b << "     " << m << "     " <<
			fx(a) << "     " << fx(b) << "     " << fx(m) << "          " << b - a << endl;
		if (fx(b) * fx(m) < 0) {
			a = m;
		}
		else {
			b = m;
		}
		k++;
		m = (a + b) / 2;
	}
	cout << k << "     " << a << "     " << b << "     " << m << "     " <<
		fx(a) << "     " << fx(b) << "     " << fx(m) << "     " << b - a << endl;

	return m;
}

double solveNewton(function fx, function dfx, double x0) {
	cout << "NEWTON" << endl;
	int k = 0;
	double x1 = x0 - fx(x0) / dfx(x0);
	cout << fixed << setprecision(12);
	cout << "k     " << "x                   " << "eps" << endl;
	while (fabs(x1 - x0) >= eps7) {
		cout << k << "     " << x0 << "     " << fabs(x1 - x0) << endl;
		x0 = x1;
		x1 = x0 - fx(x0) / dfx(x0);
		k++;
	}
	cout << k << "     " << x0 << "     " << fabs(x1 - x0) << endl;
	return x1;
}

double solveMPI(function fx, double x) {
	cout << "MPI" << endl;
	int k = 0;
	double xPrev = -fx(x);
	cout << "k     " << "x                   " << "eps" << endl;
	while (fabs(x - xPrev) >= eps7) {
		cout << k << "     " << x << "     " << fabs(x - xPrev) << endl;
		xPrev = x;
		x = fx(x);
		k++;
	}
	cout << k << "     " << x << "     " << fabs(x - xPrev) << endl;
	return x;
}

int main() {
	double d = dichotomy(fx, -5, 0);
	solveNewton(fx, dfx, d);
	solveMPI(gx, d);
	return 0;
}