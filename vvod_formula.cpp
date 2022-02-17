#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "functions.h"
using namespace std;

void vvod_formula(double*matrix, int n, int k){
    for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				matrix[i * n + j] = f(k, n, i, j ); // грамотная формула реализована
			}
		}
}
void vvod_formula2(double* matrix, int n, int k) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i * (n+1) + j] = f(k, n, i, j); // грамотная формула реализована
		}
	}
}