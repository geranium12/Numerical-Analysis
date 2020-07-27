#include "LES.h"
#include "PrintToFile.h"

void PrintToFile::operator()(LES& les, int iteration) {
	int pr, w;
	std::string nameOfFile;
	nameOfFile = les.state + ".txt";
	std::ofstream out(nameOfFile);

	if (les.state.find("RelaxationMethod") != std::string::npos) {
		pr = 6;
		w = 11;

		if (les.state.find("ExceededKMax") != std::string::npos) {
			out << "ExceededKMax\n";
		}
		else {
			out << "Iteration: " << iteration << "\n";
		};

		for (int i = 0; i < les.mSize_; ++i) {
			out << std::fixed << std::setprecision(pr) << std::setw(w)
				<< les.curSol_[i] << " \\\\\n";
		}

		out.flush();
		return;
	}

	w = 4;

	for (int i = 0; i < les.mSize_; ++i) {
		for (int j = 0; j < les.mSize_; ++j) {
			out << std::setw(w) << les[i][j] << " &";
		}

		out << "\t" << std::setw(w) << les[i][les.mSize_] << " \\\\\n";
	}

	out.flush();
}