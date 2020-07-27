#include "Matrix.h"

ExtendedMatrix::ExtendedMatrix(const ExtendedMatrix& matrix) : mSize(matrix
	.mSize) {
	a = new float* [mSize];

	for (int i = 0; i < mSize; ++i) {
		a[i] = new float[mSize + 1];

		for (int j = 0; j < mSize; ++j) {
			a[i][j] = matrix[i][j];
		}

		a[i][mSize] = matrix[i][mSize];
	}
}

ExtendedMatrix::ExtendedMatrix(int size, int solFirstVal) : mSize(size) {
	a = new float* [mSize];

	for (int i = 0; i < mSize; ++i) {
		a[i] = new float[mSize + 1];
	}

	generateCoefs();
	generateCTerms(solFirstVal);
}

ExtendedMatrix::~ExtendedMatrix() {
	for (int i = 0; i < mSize; ++i) {
		delete[] a[i];
	}
}

float* ExtendedMatrix::operator[](int index) const {
	return a[index];
}

int ExtendedMatrix::size() const {
	return mSize;
}

void ExtendedMatrix::generateCoefs() {
	std::srand(std::time(0));

	for (int i = 0; i < mSize; ++i) {
		int sum = 0;

		for (int j = 0; j < mSize; ++j) {
			a[i][j] = std::rand() % 5 - 4;
			sum += a[i][j];
		}

		sum -= a[i][i];
		a[i][i] = -sum;
	}

	a[0][0]++;
}

void ExtendedMatrix::generateCTerms(int solFirstVal) {
	for (int i = 0; i < mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			a[i][mSize] += (solFirstVal + j) * a[i][j];
		}
	}
}