#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "functions.h"
using namespace std;
int vvod_file(double*matrix, int n, string filename){ // ввод матрицы
    ifstream f(filename);
		if (!f) {
			cerr << "File could not be opened!" << endl;
			delete[] matrix;
			return -1;
		}
		for (int i = 0; i < n * n; i++) {
			if (!(f >> matrix[i])) {
				cout<<i ; 
				cerr << "Incorrect format of data!" << endl;
				delete[] matrix;
				return -1;
			}
		}
	return 0;
}