#include <iostream>
#include <vector>
#include <time.h>
#include <thread>

#include "Matrix.h"
#include "LES.h"
#include "PrintToFile.h"

const int K_MAX = 1000;
const float E = 1e-5;

int SIZE = 10;
const int RAND_RANGE = 10;
const int SETW = 10;

using namespace std;

// randomizing double number
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

// fill in Matrix A
void fillInMatrix(vector<vector<float>>& A)
{
	/*
	std::srand(std::time(0));

	for (int i = 0; i < A.size(); ++i) {
		int sum = 0;

		for (int j = 0; j < A.size(); ++j) {
			A[i][j] = std::rand() % 5 - 4;
			sum += A[i][j];
		}

		sum -= A[i][i];
		A[i][i] = -sum;
	}

	A[0][0]++;
	*/
	/*
	int value = 5;
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			A[i][A.size()] += (value + j) * A[i][j];
		}
	*/

	
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			A[i][j] = randDoubleNumber(10);
		}
	}
	
}

// fill in Vector x
void fillInVector(vector<float>& x)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] = randDoubleNumber(10);
	}
}

// input from file
void inputFromFile(istream& fin, vector<vector<float>>& A, vector<float>& f)
{
	fin >> SIZE;
	A.resize(SIZE, vector<float>(SIZE));
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			fin >> A[i][j];
		}
	}
	f.resize(SIZE);
	for (int i = 0; i < f.size(); i++)
	{
		fin >> f[i];
	}
}

// output Vector
void outputVector(ostream& fout, const vector<float>& A)
{
	for (int i = 0; i < A.size(); i++)
	{
		fout << setw(10) << setprecision(5) << A[i];
	}
	fout << endl;
}

// output Matrix
void outputMatrix(ostream& fout, const vector<vector<float>>& A, bool f = false)
{
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			if (f)
			{
				fout << fixed << setprecision(6);
			}
			fout << setw(SETW);
			fout << A[i][j];
		}
		fout << endl;
	}
}

// get result f Vector
vector<float> getResultVector(const vector<vector<float>>& A, const vector<float>& x)
{
	vector<float> f(A.size());
	for (int i = 0; i < A.size(); i++)
	{
		f[i] = 0;
		for (int j = 0; j < A.size(); j++)
		{
			f[i] += A[i][j] * x[j];
		}
		//f[i] *= 100;
		//f[i] = round(f[i]);
		//f[i] /= 100;
	}
	return f;
}

// counting max residual
float maxResidual(const vector<float>& f, const vector<float>& f1)
{
	float max = 0;
	for (int i = 0; i < f.size(); i++)
	{
		if (fabs(f[i] - f1[i]) > max)
		{
			max = fabs(f[i] - f1[i]);
		}
	}
	return max;
}

vector<float> SolRelaxationMethod(const vector<vector<float>>& A, const vector<float>& f, const double accuracy, const double option, int& k) {
	vector<float> prevSol(A.size());
	vector<float> curSol(A.size());

	int n = A.size();
	for (k = 0; k < K_MAX; ++k) {
		for (int i = 0; i < n; ++i) {
			curSol[i] = f[i];

			for (int j = 0; j < i; ++j) {
				curSol[i] -= A[i][j] * curSol[j];
			}

			for (int j = i + 1; j < n; ++j) {
				curSol[i] -= A[i][j] * prevSol[j];
			}

			curSol[i] = curSol[i] * option / A[i][i] + (1 - option) * prevSol[i];
		}

		if (maxResidual(curSol, prevSol) < accuracy) {
			return curSol;
		}

		swap(prevSol, curSol);
	}
	return vector<float>(n);
}

// multiplying Matrices
vector<vector<float>> multiplyingMatrices(const vector<vector<float>>& A1, const vector<vector<float>>& A2)
{
	vector<vector<float>> res(A1.size(), vector<float>(A1.size()));
	for (int i = 0; i < A1.size(); i++)
	{
		for (int j = 0; j < A2.size(); j++)
		{
			float temp = 0;
			for (int k = 0; k < A1[0].size(); k++)
			{
				temp += (A1[i][k] * A2[k][j]);
			}
			res[i][j] = temp;
		}
	}
	return res;
}

vector<vector<float>> transpose(const vector<vector<float>>& A)
{
	vector<vector<float>> ATranspose(A.size(), vector<float>(A.size()));
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			ATranspose[i][j] = A[j][i];
		}
	}
	return ATranspose;
}

void doSimmetr(vector<vector<float>>& A)
{
	vector<vector<float>> ATranspose = transpose(A);
	A = multiplyingMatrices(ATranspose, A);
}

int main() {
	
	// Preparing data
	vector<vector<float>> A(SIZE, vector<float>(SIZE));
	fillInMatrix(A);
	doSimmetr(A);

	vector<float> x(SIZE);
	fillInVector(x);

	vector<float> f = getResultVector(A, x);
	vector<float> f2 = f;

	vector<vector<float>> ATranspose = transpose(A);
	for (int i = 0; i < f.size(); i++)
	{
		f[i] = f[i] * A[i][0];
	}
	/*for (int i = 0; i < f.size(); i++)
	{
		f[i] = f1[0][i];
	}
	*/
	// Writing data to output.txt
	ofstream fout("output.txt");
	fout << "Matrix A:" << endl;
	outputMatrix(fout, A);

	fout << "Vector f:" << endl;
	outputVector(fout, f);

	fout << "Vector x:" << endl;
	outputVector(fout, x);
	
	double w;
	int k;
	vector<float> x1;
	f = f2;

	fout.width(10);
	fout << "w" << setw(10) << "k + 1" << setw(15) << "maxResidual" << endl;

	w = 0.2;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 0.5;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 0.8;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 1;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 1.3;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 1.5;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;

	w = 1.8;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << setprecision(2) << w;
	fout << setprecision(6) << fixed << setw(10) << k + 1 << setw(10) << maxResidual(x, x1) << endl;
	outputVector(fout, x1);
	
	/*
	fout = ofstream("input3.txt");
	fout << SIZE << endl;
	outputMatrix(fout, A);

	outputVector(fout, f);

	outputVector(fout, x);
	fout.close();

	ifstream fin("input3.txt");
	inputFromFile(fin, A, f);

	// Writing data to output.txt
	fout = ofstream("output4.txt");
	fout << "Matrix A:" << endl;
	outputMatrix(fout, A);

	fout << "Vector f:" << endl;
	outputVector(fout, f);
	fout << endl;
	
	w = 0.2;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout.width(10);
	fout << setprecision(2) << fixed;
	fout << "w" << setw(10) << "k + 1" << endl;
	w = 0.5;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;
	w = 0.8;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;

	w = 1;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;

	w = 1.3;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;

	w = 1.5;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;

	w = 1.8;
	x1 = SolRelaxationMethod(A, f, E, w, k);
	fout << w << setw(10) << k + 1 << endl;
	*/
	return 0;
}