#pragma once
#include <ctime>
#include <cstdlib>

class Matrix {
public:
	Matrix(const Matrix&);
	Matrix(int, int);
	~Matrix();

	float* operator[] (int) const;
	int size() const;

private:
	void generateCoefs();
	void generateCTerms(int);

	float** a;
	int mSize;
};