#pragma once
#include <string>
#include <cmath>

class LES {
public:

	LES(Matrix);

	int findSolRelaxationMethod(int, double, double);

private:
	Matrix matrix;
	int mSize_;
	float* curSol_;
	float* prevSol_;
	std::string state;

	friend class PrintToFile;

	float maxDiff();
	float* operator[] (int);
};