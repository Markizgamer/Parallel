#ifndef functions_H
#define FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
void vvod_formula(double* matrix, int n, int k);
void vvod_formula2(double* matrix, int n, int k);
int vvod_file(double* matrix, int n, string filename);
int jordan(double* matrix, double* b, double* x, int n);
double f(int k, int n, int i, int j);
void print_matrix(int l, int n, double* matrix, int m);
double norm(double* v, int n);
double* multiply(double* left, double* right, int n);
double diff(double* coefs, double* b, double* x, int n);
void fillB(double* coeff, double* B, int n);
void fillB2(double* coeff, double* B, int n);
#endif