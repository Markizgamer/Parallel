#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <mm_malloc.h>
#include <pthread.h>
#include <sys/time.h>


using namespace std;
void array_output(int l, int n, int m, double *a); // где мы вообще замеряем время??? 
double defeps(int n, double *a);
void pthread_barrier(int total_threads);
double get_full_time();


double defeps(int n, double *a) //погрешность через норму
{
    if(!a)
        return 0;
    double max = 0;
    double sum_str = 0;
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            sum_str = sum_str + fabs(a[i *n + j]);
        }
        if (sum_str > max) {
            max = sum_str;
        }
        sum_str = 0;
    }
    return max*10e-15;
}


int gauss(double *a,double *b,int n, int *var, double *y, int thread_num, int total_threads)
{

    int i=0,j=0,h,l,temp;
    double k;
    double max;
    double tmp;
    int kol, first_row, last_row = 0;
    double eps = defeps(n, a);
    k = 0;
    for (i=0; i<n-1; ++i)
    {
        if (thread_num == 0)
        {
            max=a[i*n+i];
            h=i;
            for (j=i; j<=n-1; ++j)
            {
                if (fabs(a[j+n*i])>fabs(max))
                {
                    h=j; //запомнили столбец максимума
                    max = a[j+n*i];
                }
            }
           /* if (fabs(max) < eps)
            {
                return -1;
            }*/
            if (h!=i) //хотим поменять местами h и i столбцы
            {
                for (j=0;j<n;++j) //в этом цикле j - строки
                {
                    k=a[h+n*j];
                    a[h+n*j]=a[i+n*j];
                    a[i+n*j]=k;
                }
                temp = var[h];  //поменяли переменные местами, запомнили
                var[h] = var[i];
                var[i] = temp;
            }
            //k=a[i+n*i]; //максимум
        }
        //cout << "\n";
        pthread_barrier(total_threads);

        k=a[i+n*i]; //максимум
        //cout << max << " " << '\n';
	if (fabs(k) < eps)
            {
                return -1;
            }

        first_row = (n-i)*thread_num; //вот это я пока не понимаю, если быть честным 
        first_row /= total_threads;
        first_row += i;
        last_row = (n-i)*(thread_num+1);
        last_row /= total_threads;
        last_row = last_row - 1 + i;

        //cout << thread_num << ' ' << i << ' ' << first_row << ' ' << last_row << '\n';

        for (j=first_row; j<=last_row; j++) // по строкам; прямой ход гаусса
        {
            if (j != i)
            {
                tmp = a[i + n * j]; //первый элемент каждой строки подматрицы
                //cout << thread_num << ' ' << i << ' ' << first_row << ' ' << last_row << ' ' << tmp << '\n';
                b[j] = b[j] - b[i] * tmp / k;

                for (l = i; l <= n - 1; l++) //по столбцам, оставить такие границы??
                {
                    a[l + n * j] = a[l + n * j] - a[l + n * i] * tmp / k;
                }
            }
        }
        //array_output(4, 4, 4, a);
        //cout << '\n';
        pthread_barrier(total_threads);

    }
    //cout << '\n';
    //array_output(1, 4, 4, b); //вектор b вывелся правильно после всех преобразований

    //array_output(1, n, 4, a);
    //обратный ход
    if (thread_num == 0)
    {
        if (fabs(a[n*n-1]) > eps)
        {
            b[n-1] = b[n-1]/a[n*n-1];
            y[n-1] = b[n-1]/a[n*n-1];
            a[n*n-1] = 1;
        } else {
            return -1;
        }
        //cout << '\n';
        //array_output(4, 4, 4, a);
        //cout << '\n';
    }
    pthread_barrier(total_threads);

    if (thread_num == 0) //что происходит с вектором b
    {
        //cout << "Вектор b до обратного хода Гаусса" << "\n";
        //array_output(1, 4, 4, b);
        for (j = n - 2; j >= 0; --j)
        {
            //cout << b[j] << "\n";
            for (i = j+1; i <= n-1; i++)
            {
                b[j] = b[j] - a[i + n * j] * b[i];
                //cout << y[j] << "\n";
            }
            b[j] = b[j] / a[j + n * j];
            //cout << b[j] << "\n";
            //cout << "Вектор b после очередной найденной переменной" << "\n";
            //array_output(1, 4, 4, b);
        }
        //cout << "Вектор b после обратного хода Гаусса" << '\n';
        //array_output(1, 4, 4, b);
    }
    pthread_barrier(total_threads);

    if (thread_num == 0)
    {
       /* for (i=0;i<n;i++)
        {
            b[i]=b[i]/a[i+i*n];
        }*/ // "вернем" решения на место
        for (int i = 0; i < n; i++)
        {
            y[i] = b[i];
        }
        for (int i = 0; i < n; i++)
        {
            b[var[i]] = y[i]; //решение хранится в b
        }
    }

    pthread_barrier(total_threads);


    //printf("\n");
    return 0;
}
