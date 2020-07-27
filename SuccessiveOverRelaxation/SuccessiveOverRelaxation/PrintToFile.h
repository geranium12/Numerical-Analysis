#pragma once
#include <string>
#include <iomanip>
#include <fstream>

class LES;

struct PrintToFile {
	void operator() (LES&, int iteration = 0);
};