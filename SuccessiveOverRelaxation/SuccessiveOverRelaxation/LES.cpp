#include "LES.h"

LES::LES(ExtendedMatrix exMatrix)
	: mSize_(exMatrix.size()), matrix(exMatrix) {

	curSol_ = new float[mSize_];
	prevSol_ = new float[mSize_];

	state = "LinearEquationsSystem";
}

int LES::findSolRelaxationMethod(int kMax, double accuracy, double option) {
	state = "RelaxationMethod" + std::to_string(int(10 * option));

	for (int k = 0; k < kMax; ++k) {
		for (int i = 0; i < mSize_; ++i) {
			curSol_[i] = matrix[i][mSize_];

			for (int j = 0; j < i; ++j) {
				curSol_[i] -= matrix[i][j] * curSol_[j];
			}

			for (int j = i + 1; j < mSize_; ++j) {
				curSol_[i] -= matrix[i][j] * prevSol_[j];
			}

			curSol_[i] = curSol_[i] * option / matrix[i][i] + (1 - option) * prevSol_[i];
		}

		std::swap(prevSol_, curSol_);

		if (maxDiff() < accuracy) {
			return k;
		}
	}

	state += "ExceededKMax";
	return -1;
}

float LES::maxDiff() {
	float mDiff = 0;

	for (int i = 0; i < mSize_; ++i) {
		if (mDiff < std::abs(prevSol_[i] - curSol_[i])) {
			mDiff = std::abs(prevSol_[i] - curSol_[i]);
		}
	}

	return mDiff;
}

float* LES::operator[] (int i) {
	return matrix[i];
}