#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "functions.h"
using namespace std;

double f(int k, int n, int i, int j) {
	if (k == 1) {

		return (n - max(i, j));
	}
	else if (k == 2) {
		return (max(i, j));
	}
	else if (k == 3) {
		return abs(i - j);
	}
	else if (k == 4) {
		return 1.0 / (i + j + 1);
	}
	return 0;
}

void print_matrix(int l, int n, double* matrix, int m) {
	int min_rows = min(m, l),
		min_cols = min(m, n);

	for (int i = 0; i < min_rows; i++) {
		for (int j = 0; j < min_cols; j++) {
			cout << fixed << setw(10) << setprecision(3) << matrix[i * n + j];
		}
		cout << endl;
	}
}


double norm(double* v, int n) {  //

	double total = 0;
	for (int i = 0; i < n; i++) {
		if (total < fabs(v[i])) {
			total = fabs(v[i]);
		}
	}

	return total;
}

double* multiply(double* left, double* right, int n) {
	double* result = new double[n];
	for (int i = 0; i < n; i++) {
		result[i] = 0;
		for (int k = 0; k < n; k++) {
			result[i] += left[i * n + k] * right[k];
		}
	}
	return result;
}

double diff(double* coefs, double* b, double* x, int n) {
	double* m = multiply(coefs, x, n); // ��� �� �������� ������� �� ���� 

	for (int i = 0; i < n; i++) {
		m[i] -= b[i];
	}

	return norm(m, n) / norm(b, n);
}

void fillB(double* coeff, double* B, int n) {

	for (int i = 0; i < n; ++i) {
		B[i] = 0;
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < (n + 1) / 2; ++j)
			B[i] += coeff[i * n + 2 * j];
	}
}

void fillB2(double* coeff, double* B, int n) {

	for (int i = 0; i < n; ++i) {
		B[i] = 0;
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < (n + 1) / 2; ++j)
			B[i] += coeff[i * (n+1) + 2 * j];
	}
}

