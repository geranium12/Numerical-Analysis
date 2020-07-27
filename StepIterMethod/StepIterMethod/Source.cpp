#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

const double EPS = 0.000001;

// output matrix
void outputMatrix(ostream& fout, const vector<vector<float>>& A, bool f = false)
{
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j <A.size(); j++)
		{
			if (f)
			{
				fout << fixed << setprecision(6);
			}
			fout << setw(12);
			fout << A[i][j];
		}
		fout << endl;
	}
}

// output Vector
void outputVector(ostream& fout, const vector<float>& A)
{
	for (int i = 0; i < A.size(); i++)
	{
		fout << A[i] << endl;
	}
}

void createMatrixA(istream& fin, vector<vector<float>>& A, int k)
{
	int n;
	fin >> n;
	A.resize(n, vector<float>(n));
	vector<vector<float>> B(n, vector<float>(n));
	for (int i = 0; i < B.size(); i++)
	{
		for (int j = 0; j < B.size(); j++)
		{
			fin >> B[i][j];
		}
	}
	vector<vector<float>> C(n, vector<float>(n));
	for (int i = 0; i < C.size(); i++)
	{
		for (int j = 0; j < C.size(); j++)
		{
			fin >> C[i][j];
		}
	}

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			A[i][j] = B[i][j] * k + C[i][j];
		}
	}
}

void inputFromFile(ifstream& fin, vector<vector<float>>& A)
{
	int n;
	fin >> n;
	A.resize(n, vector<float>(n));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			fin >> A[i][j];
		}
	}
}

// count max residual
float maxResidual(const vector<float>& f)
{
	float max = 0;
	for (int i = 0; i < f.size(); i++)
	{
		if (fabs(f[i]) > max)
		{
			max = fabs(f[i]);
		}
	}
	return max;
}

// find scalar multiplying
float scalMulti(vector<float>& y1, vector<float>& y0)
{
	float scalMulti = 0;
	for (int i = 0; i < y1.size(); i++)
	{
		scalMulti += (y1[i] * y0[i]);
	}
	return scalMulti;
}

// multiply matrices
vector<float> multiplyingMatrices(const vector<vector<float>>& A1, const vector<float>& A2)
{
	vector<float> res(A2.size());
	for (int i = 0; i < A1.size(); i++)
	{
		float temp = 0;
		for (int j = 0; j < A2.size(); j++)
		{
			temp += (A1[i][j] * A2[j]);
		}
		res[i] = temp;
	}
	return res;
}

vector<float> divVector(vector<float> f, float k)
{
	for (int i = 0; i < f.size(); i++)
	{
		f[i] /= k;
	}
	return f;
}

vector<float> divVectorCoord(vector<float> f, vector<float> f2)
{
	for (int i = 0; i < f.size(); i++)
	{
		f[i] /= f2[i];
	}
	return f;
}

vector<float> subVectorCoord(vector<float> f, vector<float> f2)
{
	for (int i = 0; i < f.size(); i++)
	{
		f[i] -= f2[i];
	}
	return f;
}

// max lambda - scalar mul
pair<float, vector<float>> findMaxLambda1(ostream& fout, vector<vector<float>>& A, vector<float> r0, int& countIter)
{
	vector<float> Ar0;
	vector<float> r1;

	float lambda0;
	float lambda1;
	float sub = 1;

	Ar0 = multiplyingMatrices(A, r0);
	r1 = divVector(Ar0, maxResidual(Ar0));
	lambda0 = scalMulti(r0, Ar0) / scalMulti(r0, r0);

	r0 = r1;
	Ar0 = multiplyingMatrices(A, r0);
	r1 = divVector(Ar0, maxResidual(Ar0));
	lambda1 = scalMulti(r0, Ar0) / scalMulti(r0, r0);
	
	sub = lambda1 - lambda0;
	countIter = 1;

	while (fabs(sub) > EPS)
	{
		lambda0 = lambda1;
		r0 = r1;
		Ar0 = multiplyingMatrices(A, r0);
		r1 = divVector(Ar0, maxResidual(Ar0));
		lambda1 = scalMulti(r0, Ar0) / scalMulti(r0, r0);
		sub = lambda1 - lambda0;
		countIter++;
	}
	return make_pair(lambda0, r0);
}

// max lambda - coord div
pair<float, vector<float>> findMaxLambda2(ostream& fout, vector<vector<float>>& A, vector<float> r0, int& countIter)
{
	vector<float> Ar0;
	vector<float> lambda0;
	vector<float> lambda1;
	vector<float> sub(A.size(), 1);

	r0 = divVector(r0, maxResidual(r0));
	Ar0 = multiplyingMatrices(A, r0);
	lambda0 = divVectorCoord(Ar0, r0);

	r0 = Ar0;
	r0 = divVector(r0, maxResidual(r0));
	Ar0 = multiplyingMatrices(A, r0);
	lambda1 = divVectorCoord(Ar0, r0);

	sub = subVectorCoord(lambda1, lambda0);
	countIter = 1;

	while (maxResidual(sub) > EPS)
	{
		lambda0 = lambda1;
		r0 = Ar0;
		r0 = divVector(r0, maxResidual(r0));
		Ar0 = multiplyingMatrices(A, r0);
		lambda1 = divVectorCoord(Ar0, r0);
		sub = subVectorCoord(lambda1, lambda0);
		countIter++;
	}

	int maxPos = 0;
	float maxLambda = lambda0[maxPos];
	for (int i = 1; i < lambda0.size(); i++)
	{
		if (fabs(lambda0[i]) > fabs(maxLambda))
		{
			maxPos = i;
			maxLambda = lambda0[maxPos];
		}
	}
	return make_pair(maxLambda, r0);
}

int main()
{
	ifstream fin("input.txt");
	ofstream fout("output.txt");
	vector<vector<float>> A;

	//inputFromFile(fin, A);
	createMatrixA(fin, A, 7);
	fout << "Matrix A : " << endl;
	outputMatrix(fout, A, true);

	fout << "1 VARIANT (Scalar mul)" << endl;
	vector<float> r0 = { 4, 3, 2, 1, 1, 1 };
	int countIter;
	pair<float, vector<float>> res = findMaxLambda1(fout, A, r0, countIter);
	
	fout << "Initial vector : " << endl;
	outputVector(fout, r0);
	fout << "Max Lambda : " << res.first << endl;
	fout << "Characteristic vector :" << endl;
	outputVector(fout, res.second);
	fout << "Number of iterations : " << countIter << endl;

	fout << "2 VARIANT (Coord div)" << endl;
	//r0.resize(A.size(), 1);
	r0 = { 4, 3, 2, 1, 1, 1 };
	res = findMaxLambda2(fout, A, r0, countIter);
	
	fout << "Initial vector : " << endl;
	outputVector(fout, r0);
	fout << "Max Lambda : " << res.first << endl;
	fout << "Characteristic vector :" << endl;
	outputVector(fout, res.second);
	fout << "Number of iterations : " << countIter << endl;
	return 0;
}