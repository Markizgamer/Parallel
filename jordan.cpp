#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "functions.h"
using namespace std;

void Print_matrix(int l, int n, double* matrix, int m) {
	int min_rows = min(m, l),
		min_cols = min(m, n);

	for (int i = 0; i < min_rows; i++) {
		for (int j = 0; j < min_cols; j++) {
			cout << fixed << setw(10) << setprecision(3) << matrix[i * n + j];
		}
		cout << endl;
	}
}
void swap_rows(double* coefs, int i1, int i2, int n, double* b) { // свапаем строчки соответственнь 
	double tmp; // вводим временную переменную 
	for (int i = 0; i < n; i++) {
		tmp = coefs[i1 * n + i]; // мы элементы одной строки меняем с элементами другой строки по сути
		coefs[i1 * n + i] = coefs[i2 * n + i];
		coefs[i2 * n + i] = tmp;
	}
	tmp = b[i1]; // то же самое делаем и с элементами правой части 
	b[i1] = b[i2]; // поменяли 
	b[i2] = tmp;
}


// ввели одномерный массив х-ов  // мы это в мейне создали т е видим по сути его  

void swap_cols(double* coefs, int j1, int j2, int n, int* xes) { //меняем колонки соотвественно 
	double tmp;
	int u;
	for (int i = 0; i < n; i++) {
		tmp = coefs[i * n + j1];
		coefs[i * n + j1] = coefs[i * n + j2]; // переставив коэффициенты столбцов надо переставить индексы иксов, которые пронумерованы от 0 до n-1 
		coefs[i * n + j2] = tmp;
	}

	u = xes[j1];
	xes[j1] = xes[j2];
	xes[j2] = u;
}
int jordan(double* matrix, double* b, double* x, int n) {
	//cout<<"первый проход"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<"\n";

	int* xes = new int[n]; // тут мы и заводим массив индексов и получаем n индексов от 0 до n-1
	for (int i = 0; i < n; i++) {
		xes[i] = i; // заполняем эти n индексов 
	}


	int i_max = 0, j_max = 0;
	double max = matrix[0];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (max < matrix[i * n + j]) {
				i_max = i;
				j_max = j;
				max = matrix[i * n + j];
			}
		}
	}
	//cout<<"максимальный элемент"<<"\n";
	//cout<<max;

	swap_rows(matrix, 0, i_max, n, b);

	swap_cols(matrix, 0, j_max, n, xes);
	//cout<<"матрица посе перестановки"<<"\n";
	//Print_matrix(n, n, matrix , n );

	//cout<<"матрица иксов посе перестановки"<<"\n";
	//for (int o=0; o<n ; o++){
	//	cout<<xes[o]<<"\n";
//	}
	//cout<<"\n";
	//cout<<"\n";
	double tmp = matrix[0];
	//cout<<"tmp= "<<tmp<<"\n";
	//cout<<"новая матрица бешек после перестановки"<<"\n";
	//Print_matrix(n, 1, b, n);
	//cout<<"\n";
	if (FP_ZERO == fpclassify(tmp)) {
		//print_matrix(n, n, matrix, n);
		cerr << "Incorrect format of data!" << endl;
		delete[] matrix;
		delete[] xes;
		delete[] b;
		return -1; // если вдруг не дай бог у нас вылезет ноль, мы прекратим цикл и не станем делить 

	}

	for (int i = 0; i < n; i++) { // сдесь мы поделим на наш начальный коэффициент 
		matrix[0 * n + i] /= tmp;
	}
	//cout<<"Наша матрица, после того как 1 строку поделили на tmp"<<"\n";
	//Print_matrix(n, n, matrix , n );

	b[0] /= tmp; // бешку тоже делим
	//cout<<"поделили на 2"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<"\n";

	double c = 0;
	for (int i = 1; i < n; i++) { // в этом цикле мы зануляем все под и преобразуем матрицу, путем вычитания первой строки домноженной на соответствующи коэффициент
		c = matrix[i * n + 0];
		for (int j = 0; j < n; j++) {
			matrix[i * n + j] -= c * matrix[0 * n + j];
		}
		//cout<<"бешка до:"<<"\n";
		//cout<<b[i]<<"\n";
		b[i] -= c * b[0];
		//cout<<"бешка после:"<<"\n";
		//cout<<b[i]<<"\n";
	}
	//cout<<"таким образом матрица бешек после первого прохода:"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<'\n';
	//cout<<"матрица посе прохода"<<"\n";
		//Print_matrix(n, n, matrix , n );
		//cout<<'\n';
	//cout<<"основной проход:"<<"\n";
	for (int k = 1; k < n; k++) {
		//cout<<"номер прохода: "<<k<<'\n';// далее цикл для остальных размерностей матриц идет у нас 
		max = matrix[k * n + k];
		i_max = k;
		j_max = k;
		for (int i = k; i < n; i++) {
			for (int j = k; j < n; j++) {
				if (max < matrix[i * n + j]) {
					i_max = i;
					j_max = j;
					max = matrix[i * n + j];
				}
			}
		}
		//cout<<"максимальный элемент"<<"\n";
		//cout<<max;

		swap_rows(matrix, k, i_max, n, b); // меняем колонки и строки и индексы там внутри 
		swap_cols(matrix, k, j_max, n, xes);
		//cout<<"матрица посе перестановки"<<"\n";
		//Print_matrix(n, n, matrix , n );
		//cout<<"\n";// передадим сюда свеже созданный массив
		//cout<<"матрица бешек посе перестановки"<<"\n";
		//Print_matrix(n, 1, b , n );
		//cout<<"\n";
		//cout<<"i_max = "<<i_max<<"j_max ="<<j_max<<"\n";
		//cout<<"матрица иксов посе перестановки"<<"\n";
	//for (int o=0; o<n ; o++){
	//	cout<<xes[o]<<"\n";
	//}
		tmp = matrix[k * n + k];
		//cout<<"tmp = "<<tmp<<"\n";

		if (FP_ZERO == fpclassify(tmp)) {
			//print_matrix(n, n, matrix, n);
			cerr << "Incorrect format of dat!" << endl;
			delete[] matrix;
			delete[] xes;
			delete[] b;
			return -1; // если вдруг не дай бог у нас вылезет ноль, мы прекратим цикл и не станем делить 

		}

		for (int i = k; i < n; i++) {
			matrix[k * n + i] /= tmp;
		}
		//cout<<"Матрица поле деления: "<<'\n'; 
		//Print_matrix(n,n,matrix,n);
		//cout<<'\n';

		b[k] /= tmp;
		//cout<<"матрица бешек после деления "<<"\n";
		//Print_matrix(n, 1, b , n );
		//cout<<'\n';
		for (int i = 0; i < n; i++) {
			if (i == k) { // тут мы делаем то же самое, но чтобы не занулить наше матрицу впринципе 
				continue;
			}
			c = matrix[i * n + k];
			//cout<<"c = "<<c<<"i = "<<i<<"\n";
			for (int j = k; j < n; j++) {
				matrix[i * n + j] -= c * matrix[k * n + j];
			}
			//cout<<"бешка до:"<<"\n";
		//cout<<b[i]<<"\n";
		//cout<<b[k]<<"\n";
			b[i] -= c * b[k];
			//cout<<"бешка после:"<<"\n";
		//cout<<b[i]<<"\n";
		}

		//cout<<'\n';
		//cout<<"Матрица: "<<'\n'; 
		//Print_matrix(n,n,matrix,n);
		//cout<<'\n';
		//cout<<"Матрица бешек: "<<'\n'; 
		//Print_matrix(n,1,b,n);
		//cout<<endl;
		//cout << endl << k << "-Matrix:" << endl; 

		 // после каждого шага будем выводить матрицу 
		//cout << endl;

	}





	//for (int i = 0; i < n; i++) {
	//	x[i] = b[i]; // после всего мы получаем иксы и заполняем нашу матрицу иксов бешками 
	//}
	/*cout<<"матрица иксов посе перестановки"<<"\n";
	for (int o=0; o<n ; o++){
		cout<<xes[o]<<"\n";
	}*/
	for (int i = 0; i < n; i++) {
		x[xes[i]] = b[i];
	}

	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (xes[j] == i) {
				x[j] = b[i];
				cout<<x[j]<<"\n";
			}
		}

	}*/
	delete[] xes;
	//	for (int k = n - 1; k >= 0; k--) {
	//		x[k] = b[k];
	//		for (int i = 0; i < k; i++) {
	//			b[i] -=  coefs[i][k] * x[k];
	//		}
	//	}
	return 0;
}