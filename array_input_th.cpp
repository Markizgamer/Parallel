#include <iostream>

using namespace std;

double f(int k, int n, int i, int j);
int array_input(int k, int n, double*a);
int array_from_file(int n, char* filename, double*a);


int array_input(int k, int n, double*a)
{
    if (f(k, n, 0, 0) == -1)
    {
        return -1;
    }
    else
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i * n + j] = f(k, n, i, j);
            }
        }
    }
    return 0;
}

double f(int k, int n, int i, int j)
{
    switch(k)
    {
        case 1:
            return n - max(i, j) + 1;
        case 2:
            return max(i, j);
        case 3:
            return max(i - j, j - i);
        case 4:
            return 1./ (i + j + 1);
        case 5:
              if (i == n/2)
                  return 0;
              else
                  return n - max(i, j) + 1;
        default:
            //printf("Такой формулы нет\n");
            return -1;
    }
}

int array_from_file(int n, char* filename, double*a)  // убрать выводы в main
{
    FILE *fin;

    fin = fopen(filename, "r");
    if (!(fin))
    {
        //printf("Файл не открылся\n");
        return -1;
    }
    for (int i = 0; i<n*n; i++)
    {
        if (((fscanf(fin, "%lf", &a[i])) != 1) && (!(feof(fin))))
        {
            //printf("Некорректные данные\n");
            fclose(fin);
            return -2;
        }
    }
    fclose(fin);
    return 0;
}


