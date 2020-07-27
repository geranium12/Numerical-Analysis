#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <chrono>
#include <thread>

using namespace std;


int SIZE = 10;
const int RAND_RANGE = 10;
const int SETW = 10;

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

// Gaussian Elimination
vector<float> gaussianElimination(vector<vector<float>> A, vector<float> f, float& det)
{
	/*
		vector<float> coeff(A.size() - j - 1);
		for (int k = 0; k < coeff.size(); k++)
		{
			coeff[k] = A[A.size() - coeff.size() + k][j] / A[j][j];
		}
		for (int k = 0; k < coeff.size(); k++)
		{
			for (int i = 0; i <= j; i++)
			{
				A[k + 1][i] = 0;
			}
			
			for (int i = 0; i < A.size(); i++)
			{
				A[k + 1][i] -= A[j][i] * coeff[k];
			}
			f[k + j + 1] -= f[j] * coeff[k];
		}
	}

	det = 0;
	for (int i = 0; i < A.size(); i++)
	{
		det *= A[i][i];
	}
	if (!count % 2)
	{
		det *= -1;
	}

	vector<float> res(A.size());
	res[A.size() - 1] = f[f.size() - 1] / A[A.size() - 1][A.size() - 1];
	for (int i = res.size() - 2; i >= 0; i--)
	{
		float temp = f[i];
		for (int j = res.size() - 1; j > i; j--)
		{
			temp -= A[i][j] * res[j];
		}
		res[i] = temp / A[i][i];
	}
	return res;
	*/

	vector<float> res(A.size());
	int n = A.size();
	int count = 0;
	for (int j = 0; j < n; j++)
	{
		// search max index
		int maxIndex = j;
		for (int k = j + 1; k < A.size(); k++)
		{
			if (fabs(A[k][j]) > fabs(A[maxIndex][j]))
			{
				maxIndex = k;
			}
		}
		
		// swap rows
		if (maxIndex != j)
		{
			swap(A[maxIndex], A[j]);
			swap(f[maxIndex], f[j]);
			count++;
		}
		
		// multiplying row to coeff
		for (int i = j + 1; i < n; i++)
		{
			float c = -A[i][j] / A[j][j];
			for (int k = j; k < n; k++)
			{
				if (j == k)
				{
					A[i][k] = 0;
				}
				else
				{
					A[i][k] += c * A[j][k];
				}
			}
			f[i] += f[j] * c;
		}
	}

	ofstream fout("output2.txt");
	fout << "Matrix A:" << endl;
	outputMatrix(fout, A);

	// determinant of Matrix A
	det = 1;
	for (int i = 0; i < n; i++)
	{
		det *= A[i][i];
	}
	if (!det)
	{
		return res;
	}
	// result Vector x1
	for (int i = n - 1; i >= 0; i--)
	{
		res[i] = f[i] / A[i][i];
		for (int k = i - 1; k >= 0; k--)
		{
			f[k] -= A[k][i] * res[i];
		}
	}
	return res;
}

// find inverse Matrix A^-1
vector<vector<float>> inverseMatrix(vector<vector<float>> A)
{
	vector<vector<float>> res(A.size(), vector<float>(A.size()));
	vector<vector<float>> f(A.size(), vector<float>(A.size()));
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			if (i == j)
			{
				f[i][j] = 1;
			}
			else
			{
				f[i][j] = 0;
			}
		}
	}

	float det;
	for (int i = 0; i < A.size(); i++)
	{
		res[i] = gaussianElimination(A, f[i], det);
	}
	 
	vector <vector<float>> ans(A.size(), vector<float>(A.size()));
	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A.size(); i++)
		{
			ans[j][i] = res[i][j];
		}
	}
	return ans;
}

// multiplying Matrices
vector<vector<float>> multiplyingMatrices(const vector<vector<float>>& A1, const vector<vector<float>>& A2)
{
	vector<vector<float>> res(A1.size(), vector<float>(A1.size()));
	for (int i = 0; i < A1.size(); i++)
	{
		for (int j = 0; j < A1.size(); j++)
		{
			float temp = 0;
			for (int k = 0; k < A1.size(); k++)
			{
				temp += (A1[i][k] * A2[k][j]);
			}
			res[i][j] = temp;
		}
	}
	return res;
}

int main()
{
	/*
	// Preparing data
	vector<vector<float>> A(SIZE, vector<float>(SIZE));
	fillInMatrix(A);

	vector<float> x(SIZE);
	fillInVector(x);

	vector<float> f = getResultVector(A, x);

	// Writing data to output.txt
	ofstream fout("output.txt");
	fout << "Matrix A:" << endl;
	outputMatrix(fout, A);

	fout << "Vector f:" << endl;
	outputVector(fout, f);

	fout << "Vector x:" << endl;
	outputVector(fout, x);

	// Gaussian elimination
	float det;
	vector<float> x1 = gaussianElimination(A, f, det);
	if (det == 0)
	{
		fout << "Gaussian elimination can't be used for this matrix!" << endl;
		return 0;
	}
	fout << "Received Vector-Solution x1:" << endl;
	outputVector(fout, x1);

	// Max Residual
	vector<float> f1 = getResultVector(A, x1);
	fout << "Max residual ||Ax1 - f|| = " << setprecision(10) << maxResidual(f, f1) << endl;
	fout << "Max residual ||x - x1|| = " << setprecision(10) << maxResidual(x, x1) << endl;
	fout << fixed << setprecision(6);

	// determinant
	fout << "Det A: " << det << endl;

	// inverse Matrix A^-1
	fout << "A^-1:" << endl;
	vector<vector<float>> invA = inverseMatrix(A);
	outputMatrix(fout, invA, true);

	// check if gaussian algorithm is right
	fout << "A * A^-1:" << endl;
	outputMatrix(fout, multiplyingMatrices(A, invA), true);
	*/
	//Input from input.txt
	
	ifstream fin("input2.txt");
	vector<vector<float>> A;
	vector<float> f;
	inputFromFile(fin, A, f);
	float det;
	vector<float> x1 = gaussianElimination(A, f, det);
	ofstream fout("output3.txt");
	fout << "Received Vector-Solution x1:" << endl;
	outputVector(fout, x1);
	fout << "DET A: " << det << endl;
	fout << "A^-1:" << endl;
	vector<vector<float>> invA = inverseMatrix(A);
	outputMatrix(fout, invA, true);
	fout << "A * A^-1:" << endl;
	outputMatrix(fout, multiplyingMatrices(A, invA), true);
	
	return 0;
}