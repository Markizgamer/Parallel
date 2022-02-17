
#include <iostream>
#include <cstdio>

using namespace std;

void array_output(int l, int n, int m, double *a);
//void output_vector(int n, int m, double* v);

void array_output(int l, int n, int m, double *a)
{
    int a1 = min(l, m); //строки
    int a2 = min(n, m); //столбцы
    for (int i = 0; i < a1; i++)
    {
        for (int j = 0; j < a2; j++)
        {
            printf(" %10.3e", a[i*n+j]);
        }
        printf("\n");
    }
}
