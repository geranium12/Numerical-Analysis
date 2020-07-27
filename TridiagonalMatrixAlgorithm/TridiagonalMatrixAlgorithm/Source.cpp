#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <iomanip>
#include <chrono>
#include <thread>
#include <cmath>

using namespace std;

// size of generated matrix
int SIZE = 10;
// range for rand, for instance, for 10 it will be like [-10, 10]
const int RAND_RANGE = 10;
const int SETW = 10;

// function which returns random double number
float randDoubleNumber(const int a, const int b = 0)
{
	srand(time(NULL));
	std::chrono::seconds dura(1);
	std::this_thread::sleep_for(dura);
	float d = float(rand() % (a - 1)) + b;
	d *= 100;
	d += rand() % 100;
	d /= 100;
	if ((rand() % 25) % 2)
	{
		d *= -1;
	}
	return d;
}

// fill in vectors A, B, C
void fillInMatrixRandom(vector<float>& A, vector<float>& B, vector<float>& C)
{
	for (int i = 0; i < A.size(); i++)
	{
		A[i] = randDoubleNumber(RAND_RANGE) * -1;

		B[i] = randDoubleNumber(RAND_RANGE) * -1;
	}
	A[0] = 0;
	B[B.size() - 1] = 0;
	C[0] = randDoubleNumber(RAND_RANGE, round(fabs(B[0])) + 1);
	for (int i = 1; i < C.size() - 1; i++)
	{
		C[i] = randDoubleNumber(RAND_RANGE, round(fabs(A[i]) + fabs(B[i])) + 2);
	}
	C[C.size() - 1] = randDoubleNumber(RAND_RANGE, round(fabs(A[A.size() - 1])) + 1);
}

// fill in accurate solution x
void fillInVectorRandom(vector<float>& x)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] = randDoubleNumber(RAND_RANGE);
	}
}

// input from file for testing this program
void inputFromFile(istream& fin, vector<float>& A, vector<float>& B, vector<float>& C, vector<float>& f, vector<float>& x)
{
	fin >> SIZE;
	A.resize(SIZE);
	for (int i = 0; i < A.size(); i++)
	{
		fin >> A[i];
	}
	C.resize(SIZE);
	for (int i = 0; i < C.size(); i++)
	{
		fin >> C[i];
		C[i] *= -1;
	}
	B.resize(SIZE);
	for (int i = 0; i < B.size(); i++)
	{
		fin >> B[i];
	}
	f.resize(SIZE);
	for (int i = 0; i < f.size(); i++)
	{
		fin >> f[i];
	}
	x.resize(SIZE);
	for (int i = 0; i < x.size(); i++)
	{
		fin >> x[i];
	}
}

// output Input Matrix
void outputInputMatrix(ostream& fout, vector<float> A, vector<float> B, vector<float> C)
{
	fout << setw(SETW);
	fout << C[0];
	fout << setw(SETW);
	fout << B[0];
	for (int i = 2; i < SIZE; i++)
	{
		fout << setw(SETW);
		fout << 0;
	}
	fout << endl;
	for (int i = 1; i < SIZE - 1; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			if (i == j)
			{
				fout << setw(SETW);
				fout << C[i];
			}
			else if (i - j == 1)
			{
				fout << setw(SETW);
				fout << A[i];
			}
			else if (j - i == 1)
			{
				fout << setw(SETW);
				fout << B[i];
			}
			else
			{
				fout << setw(SETW);
				fout << 0;
			}
		}
		fout << endl;
	}

	for (int i = 0; i < SIZE - 2; i++)
	{
		fout << setw(SETW);
		fout << 0;
	}
	fout << setw(SETW);
	fout << A[SIZE - 1];
	fout << setw(SETW);
	fout << C[SIZE - 1];
	fout << endl;
}

// output Alpha-Beta Matrix
void outputAlphaBetaMatrix(ostream& fout, vector<float> alpha, vector<float> beta)
{
	fout << setprecision(6);
	fout << setw(SETW);
	fout << 1;
	fout << setw(SETW);
	fout << alpha[1];
	for (int i = 2; i < SIZE; i++)
	{
		fout << setw(SETW);
		fout << 0;
	}
	fout << setw(SETW);
	fout << beta[1];
	fout << endl;
	for (int i = 2; i < SIZE; i++)
	{
		for (int j = 1; j < SIZE + 1; j++)
		{
			if (i == j)
			{
				fout << setw(SETW);
				fout << 1;
			}
			else if (j - i == 1)
			{
				fout << setw(SETW);
				fout << -alpha[i];
			}
			else
			{
				fout << setw(SETW);
				fout << 0;
			}
		}
		fout << setw(SETW);
		fout << beta[i];
		fout << endl;
	}

	for (int i = 0; i < SIZE - 1; i++)
	{
		fout << setw(SETW);
		fout << 0;
	}
	fout << setw(SETW);
	fout << 1;
	fout << setw(SETW);
	fout << beta[SIZE];
	fout << endl;
}

// output vector
void outputVector(ostream& fout, vector<float> A, bool f = false)
{
	if (f)
	{
		fout << setprecision(6);
	}
	for (int i = 0; i < A.size(); i++)
	{
		fout << setw(SETW);
		fout << A[i];
	}
	fout << endl;
}

// get result vector f by multiplying
vector<float> getResultVector(vector<float>& A, vector<float>& B, vector<float>& C, vector<float>& x)
{
	int size = x.size();
	vector<float> f(size);
	f[0] = -C[0] * x[0] + B[0] * x[1];
	for (int i = 1; i < size - 1; i++)
	{
		f[i] = A[i] * x[i - 1] - C[i] * x[i] + B[i] * x[i + 1];
	}
	f[size - 1] = A[size - 1] * x[size - 2] - C[size - 1] * x[size - 1];
	return f;
}

// Tridiagonal Matrix Algorithm using formulas from the report
vector<float> tridiagonalMatrixAlgorithm(ostream& fout, vector<float> A, vector<float> B, vector<float> C, vector<float> f, float& det)
{
	int size = f.size();
	vector<float> res(size);
	vector<float> alpha(size);
	vector<float> beta(size + 1);
	alpha[1] = B[0] / C[0];
	beta[1] = -f[0] / C[0];
	det = C[0];
	for (int i = 1; i < size - 1; i++)
	{
		det *= C[i] - A[i] * alpha[i];
		alpha[i + 1] = B[i] / (C[i] - A[i] * alpha[i]);
		beta[i + 1] = (A[i] * beta[i] - f[i]) / (C[i] - A[i] * alpha[i]);
	}
	det *= C[size - 1] - A[size - 1] * alpha[size - 1];
	beta[size] = (A[size - 1] * beta[size - 1] - f[size - 1]) / (C[size - 1] - A[size - 1] * alpha[size - 1]);
	
	fout << "Alpha-Beta Matrix:" << endl;
	outputAlphaBetaMatrix(fout, alpha, beta);
	res[size - 1] = beta[size];
	for (int i = SIZE - 2; i >= 0; i--)
	{
		res[i] = alpha[i + 1] * res[i + 1] + beta[i + 1];
	}
	return res;
}

// max residual
double maxResidual(vector<float> f, vector<float> f1)
{
	double max = 0;
	for (int i = 0; i < f.size(); i++)
	{
		if (fabs(f[i] - f1[i]) > max)
		{
			max = fabs(f[i] - f1[i]);
		}
	}
	return max;
}

// check if algorithm can be used or not
bool algorithmCanBeUsed(const vector<float>& A, const vector<float>& B, const vector<float>& C)
{
	if (fabs(C[0]) <= 0 || fabs(C[C.size() - 1]) <= 0 || fabs(B[0]) <= 0 ||
		fabs(A[A.size() - 1]) <= 0 || 
		(fabs(C[0]) < fabs(B[0])) || 
		(fabs(C[C.size() - 1]) < fabs(A[A.size() - 1])))
	{
		return false;
	}
	for (int i = 1; i < A.size() - 1; i++)
	{
		if (fabs(A[i]) <= 0 || fabs(B[i]) <= 0 || (fabs(C[i]) < (fabs(A[i]) + fabs(B[i]))))
		{
			return false;
		}
	}
	float alpha = B[0] / C[0];
	for (int i = 1; i < C.size() - 1; i++)
	{
		if (!(C[i] - A[i] * alpha))
		{
			return false;
		}
		alpha = B[i] / (C[i] - A[i] * alpha);
	}
	return true;
}

int main()
{
	vector<float> A(SIZE);
	vector<float> B(SIZE);
	vector<float> C(SIZE);
	fillInMatrixRandom(A, B, C);
	vector<float> x(SIZE);
	fillInVectorRandom(x);
	vector<float> f = getResultVector(A, B, C, x);

	ofstream fout("output.txt");
	fout << "Matrix A:" << endl;
	outputInputMatrix(fout, A, B, C);
	fout << "Vector x:" << endl;
	outputVector(fout, x);
	fout << "Vector f:" << endl;
	outputVector(fout, f);

	if (!algorithmCanBeUsed(A, B, C))
	{
		fout << "The Tridiagonal Matrix Algorithm can't be used for such input data!" << endl;
		return 0;
	}

	float det;
	vector<float> x1 = tridiagonalMatrixAlgorithm(fout, A, B, C, f, det);
	fout << "Received Vector x1:" << endl;
	outputVector(fout, x1, true);
	vector<float> f1 = getResultVector(A, B, C, x1);
	fout << setprecision(10);
	fout << "Max residual ||Ax1 - f|| = " << maxResidual(f, f1) << endl;
	fout << "Max residual ||x - x1|| = " << maxResidual(x, x1) << endl;
	fout << "DET A: " << det << endl;
	
	// From input.txt
	/*
	ifstream fin("input.txt");
	vector<float> A;
	vector<float> B;
	vector<float> C;
	vector<float> f;
	vector<float> x;
	inputFromFile(fin, A, B, C, f, x);
	if (!algorithmCanBeUsed(A, B, C))
	{
		cout << "The Tridiagonal Matrix Algorithm can't be used for such input data!" << endl;
		return 0;
	}

	float det;
	vector<float> x1 = tridiagonalMatrixAlgorithm(cout, A, B, C, f, det);
	cout << "Received Vector x1:" << endl;
	outputVector(cout, x1, true);
	vector<float> f1 = getResultVector(A, B, C, x1);
	cout << setprecision(10);
	cout << "Max residual ||Ax1 - f|| = " << maxResidual(f, f1) << endl;
	cout << "Max residual ||x - x1|| = " << maxResidual(x, x1) << endl;
	cout << "DET A: " << det << endl;
	*/
	return 0;
}