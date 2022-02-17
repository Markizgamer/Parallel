#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <mm_malloc.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>

//осталось переписать самого гаусса для тредов

using namespace std;

int array_input(int k, int n, double*a);
int array_from_file(int n, char* filename, double*a);
void array_output(int l, int n, int m, double *a);
int gauss(double *a,double *b,int n, int *var, double *y, int thread_num, int total_threads);
void *gauss_threated(void*pa);
double euclid(int n, double *v);
double norma_nev(int n, double *a, double *b, double *x, double* y);
double norma_pogr(int n, double *b, double* y);
double defeps(int n, double *a);
double get_full_time();
double get_time(); //для одного процесса

typedef struct _ARGS
{
    double *a; //сама матрица
    double *b; //вектор свободных членов
    int n; // размер матрицы и векторов
    int *var; //вектор для возвращения переменных на свои позиции
    double *y; // вспомогательный вектор для возвращения переменных на свои позиции
    int thread_num; // номер потока 
    int total_threads; // всего потоков 
    double time;
    int err;
} ARGS;

void *gauss_threated(void *pa) {
    ARGS *pargs = (ARGS*)pa;
    long int t;
    int i;
    printf("Thread %d started\n", pargs->thread_num), '\n';
    pthread_barrier(pargs->total_threads);
    t = get_full_time();
    pargs->err = gauss(pargs->a, pargs->b, pargs->n, pargs->var, pargs->y, pargs->thread_num, pargs->total_threads);
    pthread_barrier(pargs->total_threads);
    //printf("Thread %d finished, time = %lf\n", pargs->thread_num, t, '\n');
    pargs->time = get_full_time() - t;
    return 0;
}


int main(int argc, char* argv[])
{
    double *a, *b, *y, *x, *a1; // в b потом решение (после гаусса), в х - свободные члены, y - вспомогательный
    int *var; // чтобы поменять местами переменные в конце
    int n, m, k, total_threads;
    pthread_t *threads;
    ARGS *args; // массив для созданных потоков 
    double t_full; // время работы полное 
    clock_t start, end;


    if ((argc != 5) && (argc != 6)) //параметров стало на один больше
    {
        printf("Некорректное число параметров\n");
        return -1;
    }

    n = atoi(argv[1]);
    m = atoi(argv[2]);
    k = atoi(argv[3]);
    total_threads = atoi(argv[4]); //теперь название файла на последнем месте

    /*a = (double*)malloc(n*n*sizeof(double));
    a1 = (double*)malloc(n*n*sizeof(double)); //чтобы в норме исходная матрица
    b = (double*)malloc(n*sizeof(double));
    x = (double*)malloc(n*sizeof(double));
    y = (double*)malloc(n*sizeof(double));
    var = (int*)malloc(n*sizeof(int));
    threads = (pthread_t*)malloc(sizeof(pthread_t)*total_threads);
    args = (ARGS*)malloc(sizeof(ARGS)*total_threads);*/
    if (!(a = (double*)malloc(n*n*sizeof(double))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -3;
    }
    if (!(a1 = (double*)malloc(n*n*sizeof(double))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -4;
    }
    if (!(b = (double*)malloc(n*sizeof(double))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -5;
    }

    if (!(x = (double*)malloc(n*sizeof(double))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -6;
    }

    if (!(y = (double*)malloc(n*sizeof(double))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -7;
    } // здесь создаем массивы данных соответственно на каждом потоке 

    if (!(var = (int*)malloc(n*sizeof(int))))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -8;
    }

    if (!(threads = (pthread_t*)malloc(sizeof(pthread_t)*total_threads)))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -9;
    }

    if (!(args = (ARGS*)malloc(sizeof(ARGS)*total_threads)))
    {
        cout << "Недостаточно памяти" << '\n';
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -10;
    }



    if (((k == 0) && (argc != 6)) || ((k != 0 ) && (argc != 5))) //параметров стало на один больше
    {
        printf("Некорректные параметры\n");
        free(a);
        free(a1);
        free(b);
        free(x);
        free(y);
        free(var);
        free(threads);
        free(args);
        return -11;
    }

    for (int i = 0; i < n; i++)
    {
        var[i] = i;// что такое var 
    }

    if (k!=0) // заполнили массив
    {
        array_input(k, n, a);
        if (array_input(k, n, a) == -1)
        {
            printf("Такой формулы нет\n");
            free(a);
            free(a1);
            free(b);
            free(x);
            free(y);
            free(var);
            free(threads);
            free(args);
            return -12;
        }
        array_input(k, n, a1);
    } else {
        array_from_file(n, argv[5], a);
        switch (array_from_file(n, argv[5], a))
        {
            case 0:
                break;
            case -1:
                printf("Файл не открылся\n");
                free(a);
                free(a1);
                free(b);
                free(x);
                free(y);
                free(var);
                free(threads);
                free(args);
                return -1;
            case -2:
                printf("Некорректные данные\n");
                free(a);
                free(a1);
                free(b);
                free(x);
                free(y);
                free(var);
                free(threads);
                free(args);
                return -2;
        }
        array_from_file(n, argv[4], a1);
    }
    // заполним столбец b
    for (int i = 0; i < n; i++)
    {
        b[i] = 0;
        x[i] = 0;
        for (int j = 0; j < (n + 1) / 2; j++)
        {
            b[i] = b[i] + a[i * n + (2 * j)];
            x[i] = x[i] + a[i * n + (2 * j)];
        }
    }
    // далее выведем массив a на экран
    array_output(n, n, m, a);
    printf("\n");
    array_output(1, n, m, b); // вывели вектор
    //printf("\n");

    //создадим все для тредов

    for (int i = 0; i < total_threads; i++) // инициализация аргументов задачи ПЕРЕПИСАТЬ
    {
        args[i].a = a;
        args[i].b = x;
        args[i].n = n;
        args[i].var = var;// здесь инициализация аргументов происходит 
        args[i].y = y;
        args[i].thread_num = i;
        args[i].total_threads = total_threads;
    } // вот тут происходит процедура инициализации тредов 

    for (int i = 0; i < total_threads; i++)
    {
        if(pthread_create(threads + i, 0, gauss_threated, args + i))
        {
            cout << "Невозможно создать задачи" << '\n';
            free(a);
            free(a1);
            free(b);
            free(x);
            free(y);
            free(var);
            free(threads);
            free(args);
            return -13;
        }
    }


    for (int i = 0; i < total_threads; i++)
    {
        if(pthread_join(threads[i], 0)) // что значит вот это? 
	{
            cout << "Невозможно синхронизировать" << '\n';
            free(a);
            free(a1);
            free(b);
            free(x);
            free(y);
            free(var);
            free(threads);
            free(args);
            return -14;
        }
    }
	
    for (int i = 0; i < total_threads; i++)
    {
	if(args[i].err==-1)// если где-то у нас будет ошибка
        {
            cout << "Матрица вырожденная" << '\n';
            free(a);
            free(a1);
            free(b);
            free(x);
            free(y);
            free(var);
            free(threads);
            free(args);
            return -14;
        }
    }

    //cout << '\n';
    //array_output(n, n, m, a); //отладочный
    //cout << '\n';
    array_output(1, n, m, args[0].b);
    printf("Execution time = %lf\n", args[0].time);
    cout << '\n';
    printf("Norma pogr = %10.3e \n", norma_pogr(n, args[0].b, args[0].y));
    cout << '\n';
    //array_output(1, n, m, b);
    printf("Norma nevyazki = %10.3e \n", norma_nev(n, a1, x, b, y));
    free(a);
    free(a1);
    free(b);
    free(x);
    free(y);
    free(var);
    free(threads);
    free(args);
    return 0;
}

double euclid(int n, double *v)
{
    double a = 0;
    for (int i = 0; i < n; i++)
    {
        a = a + v[i] * v[i];
    }
    return sqrt(a);
}

double norma_nev(int n, double *a, double *b, double *x, double* y)
{
    double norma;
    for(int i = 0; i < n; i++)
    {
        y[i] = 0;
        for(int j = 0; j < n; j++)
        {
            y[i] = y[i] + a[i * n + j] * b[j];
        }
        y[i] = y[i] - x[i];
    }
    norma = euclid(n, y) / euclid(n, x);
    return norma;
}

double norma_pogr(int n, double *b, double* y)
{
    double norma;
    for(int i = 0; i < n; i++)
        y[i] = b[i] - (i + 1) % 2;
    norma = euclid(n, y);
    return norma;
}

double get_full_time() {
    struct timeval st;
    gettimeofday(&st, NULL);
    return st.tv_sec + (st.tv_usec) / 1000000.0;
}

double get_time()
{
    struct rusage buf;
    getrusage(RUSAGE_SELF, &buf);
    return buf.ru_utime.tv_sec + (buf.ru_utime.tv_usec) / 1000000.0;
}
