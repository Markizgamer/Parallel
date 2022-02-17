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
void swap_rows(double* coefs, int i1, int i2, int n, double* b) { // ������� ������� �������������� 
	double tmp; // ������ ��������� ���������� 
	for (int i = 0; i < n; i++) {
		tmp = coefs[i1 * n + i]; // �� �������� ����� ������ ������ � ���������� ������ ������ �� ����
		coefs[i1 * n + i] = coefs[i2 * n + i];
		coefs[i2 * n + i] = tmp;
	}
	tmp = b[i1]; // �� �� ����� ������ � � ���������� ������ ����� 
	b[i1] = b[i2]; // �������� 
	b[i2] = tmp;
}


// ����� ���������� ������ �-��  // �� ��� � ����� ������� � � ����� �� ���� ���  

void swap_cols(double* coefs, int j1, int j2, int n, int* xes) { //������ ������� ������������� 
	double tmp;
	int u;
	for (int i = 0; i < n; i++) {
		tmp = coefs[i * n + j1];
		coefs[i * n + j1] = coefs[i * n + j2]; // ���������� ������������ �������� ���� ����������� ������� �����, ������� ������������� �� 0 �� n-1 
		coefs[i * n + j2] = tmp;
	}

	u = xes[j1];
	xes[j1] = xes[j2];
	xes[j2] = u;
}
int jordan(double* matrix, double* b, double* x, int n) {
	//cout<<"������ ������"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<"\n";

	int* xes = new int[n]; // ��� �� � ������� ������ �������� � �������� n �������� �� 0 �� n-1
	for (int i = 0; i < n; i++) {
		xes[i] = i; // ��������� ��� n �������� 
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
	//cout<<"������������ �������"<<"\n";
	//cout<<max;

	swap_rows(matrix, 0, i_max, n, b);

	swap_cols(matrix, 0, j_max, n, xes);
	//cout<<"������� ���� ������������"<<"\n";
	//Print_matrix(n, n, matrix , n );

	//cout<<"������� ����� ���� ������������"<<"\n";
	//for (int o=0; o<n ; o++){
	//	cout<<xes[o]<<"\n";
//	}
	//cout<<"\n";
	//cout<<"\n";
	double tmp = matrix[0];
	//cout<<"tmp= "<<tmp<<"\n";
	//cout<<"����� ������� ����� ����� ������������"<<"\n";
	//Print_matrix(n, 1, b, n);
	//cout<<"\n";
	if (FP_ZERO == fpclassify(tmp)) {
		//print_matrix(n, n, matrix, n);
		cerr << "Incorrect format of data!" << endl;
		delete[] matrix;
		delete[] xes;
		delete[] b;
		return -1; // ���� ����� �� ��� ��� � ��� ������� ����, �� ��������� ���� � �� ������ ������ 

	}

	for (int i = 0; i < n; i++) { // ����� �� ������� �� ��� ��������� ����������� 
		matrix[0 * n + i] /= tmp;
	}
	//cout<<"���� �������, ����� ���� ��� 1 ������ �������� �� tmp"<<"\n";
	//Print_matrix(n, n, matrix , n );

	b[0] /= tmp; // ����� ���� �����
	//cout<<"�������� �� 2"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<"\n";

	double c = 0;
	for (int i = 1; i < n; i++) { // � ���� ����� �� �������� ��� ��� � ����������� �������, ����� ��������� ������ ������ ����������� �� �������������� �����������
		c = matrix[i * n + 0];
		for (int j = 0; j < n; j++) {
			matrix[i * n + j] -= c * matrix[0 * n + j];
		}
		//cout<<"����� ��:"<<"\n";
		//cout<<b[i]<<"\n";
		b[i] -= c * b[0];
		//cout<<"����� �����:"<<"\n";
		//cout<<b[i]<<"\n";
	}
	//cout<<"����� ������� ������� ����� ����� ������� �������:"<<"\n";
	//Print_matrix(n, 1, b , n );
	//cout<<'\n';
	//cout<<"������� ���� �������"<<"\n";
		//Print_matrix(n, n, matrix , n );
		//cout<<'\n';
	//cout<<"�������� ������:"<<"\n";
	for (int k = 1; k < n; k++) {
		//cout<<"����� �������: "<<k<<'\n';// ����� ���� ��� ��������� ������������ ������ ���� � ��� 
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
		//cout<<"������������ �������"<<"\n";
		//cout<<max;

		swap_rows(matrix, k, i_max, n, b); // ������ ������� � ������ � ������� ��� ������ 
		swap_cols(matrix, k, j_max, n, xes);
		//cout<<"������� ���� ������������"<<"\n";
		//Print_matrix(n, n, matrix , n );
		//cout<<"\n";// ��������� ���� ����� ��������� ������
		//cout<<"������� ����� ���� ������������"<<"\n";
		//Print_matrix(n, 1, b , n );
		//cout<<"\n";
		//cout<<"i_max = "<<i_max<<"j_max ="<<j_max<<"\n";
		//cout<<"������� ����� ���� ������������"<<"\n";
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
			return -1; // ���� ����� �� ��� ��� � ��� ������� ����, �� ��������� ���� � �� ������ ������ 

		}

		for (int i = k; i < n; i++) {
			matrix[k * n + i] /= tmp;
		}
		//cout<<"������� ���� �������: "<<'\n'; 
		//Print_matrix(n,n,matrix,n);
		//cout<<'\n';

		b[k] /= tmp;
		//cout<<"������� ����� ����� ������� "<<"\n";
		//Print_matrix(n, 1, b , n );
		//cout<<'\n';
		for (int i = 0; i < n; i++) {
			if (i == k) { // ��� �� ������ �� �� �����, �� ����� �� �������� ���� ������� ��������� 
				continue;
			}
			c = matrix[i * n + k];
			//cout<<"c = "<<c<<"i = "<<i<<"\n";
			for (int j = k; j < n; j++) {
				matrix[i * n + j] -= c * matrix[k * n + j];
			}
			//cout<<"����� ��:"<<"\n";
		//cout<<b[i]<<"\n";
		//cout<<b[k]<<"\n";
			b[i] -= c * b[k];
			//cout<<"����� �����:"<<"\n";
		//cout<<b[i]<<"\n";
		}

		//cout<<'\n';
		//cout<<"�������: "<<'\n'; 
		//Print_matrix(n,n,matrix,n);
		//cout<<'\n';
		//cout<<"������� �����: "<<'\n'; 
		//Print_matrix(n,1,b,n);
		//cout<<endl;
		//cout << endl << k << "-Matrix:" << endl; 

		 // ����� ������� ���� ����� �������� ������� 
		//cout << endl;

	}





	//for (int i = 0; i < n; i++) {
	//	x[i] = b[i]; // ����� ����� �� �������� ���� � ��������� ���� ������� ����� ������� 
	//}
	/*cout<<"������� ����� ���� ������������"<<"\n";
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