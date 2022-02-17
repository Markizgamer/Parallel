#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cstdio>
#include <cmath>
#include "mpi.h"
#include "functions.h"

using namespace std;


int main(int argc, char** argv) {
	//srand(time(nullptr));
	int rank = 0;
	int nprocs = 0;
	MPI_Init(&argc, &argv);


	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n, m, k;
	if (rank == 0) {

		std::cout << "nproc = " << nprocs << " rank = " << rank << std::endl;
		
		for (int i = 0; i < argc; ++i) {
			std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
		}
		
		
		if (argc < 4 || argc > 5) {  // 
			cerr << "Invalid count of arguments" << endl;
			return 0;
		}
		

		n = stoi(argv[1]); //
		m = stoi(argv[2]);
		k = stoi(argv[3]);



		if (nprocs > n) {

			double* A; // 
			double* A2;
			A = new double[n * n];
			A2 = new double[n * n];

			string filename;
			int a;
			if (argc == 5) {

				filename = argv[4];
				a = vvod_file(A, n, filename);

				if (a != 0) {
					return -1;
				}
			}
			else {
				std::cout << "vvod formula" << std::endl;
				vvod_formula(A, n, k);
			}

			for (int j = 0; j < n * n; ++j) {
				A2[j] = A[j];
			}

			double* B = new double[n];
			fillB(A, B, n);

			double* X = new double[n];



			print_matrix(m, n, A, m);

			cout << endl;

			cout << "B: " << endl;
			for (int i = 0; i < n; i++) {
				cout << fixed << setw(10) << setprecision(3) << B[i];
			}
			cout << endl;



			clock_t begin = clock();
			a = jordan(A, B, X, n);
			if (a != 0) {
				return -1;
			}
			clock_t end = clock();

			double t = 1.0 * (end - begin) / CLOCKS_PER_SEC;

			cout << "Time: " << t << " sec." << endl;
			cout << endl;
			cout << X[0] << "\n";
			//print_matrix(n, n , coefs, m );

			cout << "Solution: " << endl;

			for (int i = 0; i < n; i++) {
				cout << fixed << setw(10) << setprecision(3) << X[i];
			}
		}
		else {

			double* AB = new double[n * (n + 1)]; // если у нас число тредов меньше n то есть смысл запускать их 
			vvod_formula2(AB, n, k);
			double* B = new double[n];
			double* R = new double[n + 1];
			fillB2(AB, B, n);
			for (int i = 0; i < n; ++i) {
				AB[i * (n + 1) + n] = B[i];
			}
			//char bufname[16] = {0};
			//snprintf(bufname, 16, "result.txt");

			//FILE* pf = fopen(bufname, "w");

			

		/*	for (int row = 0; row < n; ++row) { //тут происходит заполнение файла?
				for (int col = 0; col <= n; ++col) {
					fprintf(pf, "%.6lf  ", AB[row * (n + 1) + col]);
				}
				fprintf(pf, "\n");
			}
			fprintf(pf, "\n-------------------------------\n");*/

			time_t tStart = clock();// это необходимо переписать 

			// first step
			int bestRow = 0; // индекс лучшей строки первый щаг сделаем прямо в главной нити 
			for (int i = 1; i < n; i++) {
				if (fabs(AB[i * (n + 1)] > fabs(AB[bestRow * (n + 1)]))) {
					bestRow = i;
				}
			}
			//fprintf(stdout, "bestrow = %d\n", bestRow); // будет прям выводть лучшую строку 
			if (bestRow != 0) {
				for (int j = 0; j < n + 1; ++j) {
					double tmp = AB[j];
					AB[j] = AB[bestRow * (n + 1) + j];
					AB[bestRow * (n + 1) + j] = tmp;
				}
			}

			double A11 = AB[0]; // вот тут мы будем первую строчку матрицы препарировать 
			for (int j = 0; j < n + 1; ++j) {

				AB[j] /= A11;
			}

			double* X = new double[n];

			//

			int nrows = n / nprocs; //   здесь создаем число передаваемых строк для элементов из побочных нитей
			int rootrows = n - (nprocs - 1) * nrows; // здесь получаются оставшиеся нити  как я понял мы все равно первые строки передаем нашец мастер нити 
					
			
			for (int i = 1; i < nprocs; ++i) {
				int offset = (n+1) * (rootrows + (i-1) * nrows); // задает смещение, то есть offset- это смещение, AB изначально будет указывать на начло массива
			//	std::cout << "offset = " << offset << std::endl;
				MPI_Send(AB + offset, nrows * (n+1), MPI_DOUBLE, i, 0, MPI_COMM_WORLD); // передаем это всем тредам за исключением нашего главного
			}


			
			for (int i = 1; i < nprocs; ++i) {
				MPI_Send(AB, (n + 1), MPI_DOUBLE, i, 0, MPI_COMM_WORLD); // передаем первую строчку матрицы всем 
			}


			
			/*for (int row = 0; row < rootrows; ++row) {  // 
				fprintf(pf, "A[%d] = ", row);  // записываем в файл
				for (int col = 0; col < n+1; ++col) {
					fprintf(pf, "%.5lf  ", AB[row * (n+1) + col]);
				}
				fprintf(pf, "\n"); // здесь идет некоторая запись в файл 
			}*/

			// 
			for (int row = 1; row < rootrows; ++row) { // 
				double Ai1 = AB[row * (n + 1)];// i я строка 1 столбец фиксированный 
				for (int j = 0; j < n + 1; ++j) { // 
					AB[row * (n + 1) + j] -= Ai1 * AB[j]; // Ab[j] - это первая строчка нашей матрицы мы домножаем ее на соответству
				}
			}

			
			double* bestElemPerThread = new double[nprocs]; // это все еще мы аходимся в мастер нити  заводитм переменную в которой буут содержаться лучшие эементы найденные в каждой нити 

			//  k - это число проходов по матрице и выполнения там алгоритмических дейфствий 
			for (int k = 1; k < n; ++k) {
				//fprintf(pf, "----------------------------------------\n");
				//std::cout << "k = " << k << std::endl;

				// 
				int currentThread = 0; // текущий номер треда 

				if (k >= rootrows) { // если k будет больше числа строк, которые у нас в мастере лежат, что изменяем 
					currentThread = 1 + (k - rootrows) / nrows; // то изменяем 
				}
				
				//fprintf(pf, "\nk=%d----------------------------------------\n", k);
				//fprintf(pf, "current thread = %d\n", currentThread);
				//fflush(pf);
				/*
				for (int row = 0; row < rootrows; ++row) {
					fprintf(pf, "AB[%d] = ",  row);
					for (int col = 0; col < n + 1; ++col) {
						fprintf(pf, "%.5lf  ", AB[row * (n + 1) + col]);
					}
					fprintf(pf, "\n");
				}
				fflush(pf);
				*/

				// find aik with max abs 
				int curBestRow = 0; // текущая лучшая строка индекс
				bestElemPerThread[0] = 0;  // матрица элементов 
				if (k < rootrows) { // если k включает мастертред то ищем прм в нем наибольший элемент 
					for (int i = k; i < rootrows; ++i) {

						if (fabs(AB[i * (n + 1) + k]) > fabs(bestElemPerThread[0])) {
							bestElemPerThread[0] = AB[i * (n + 1) + k];
							curBestRow = i;
						}
					}
				}
				// receive best from each thread
				for (int th = 1; th < nprocs; ++th) {
					MPI_Recv(bestElemPerThread + th, 1, MPI_DOUBLE, th, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE); // число объекор забиваем наш буфер 
				} // вот таким образом мы забиваем наш буффер 
				//std::cout << "Received best elems" << std::endl;
				/*
				for (int i = 0; i < nprocs; ++i) {
					std::cout << bestElemPerThread[i] << std::endl;
				}
				*/
				int thWithMaxElem = currentThread; // пускай у нас максимальный элемент сначала в мастер треде 
				for (int i = currentThread; i < nprocs; ++i) {
					if (fabs(bestElemPerThread[i]) > fabs(bestElemPerThread[thWithMaxElem])) {
						thWithMaxElem = i;  // находим победителя, индекс треда с максимальным
					}
				}


				if (fabs(bestElemPerThread[thWithMaxElem]) < 1.0E-15) { // если победитель нулевой то тогда выводим ошибку по хорошему нам память надо чистить 

					std::cout << "Wrong data!!!" << std::endl;
					MPI_Finalize();
					return 1;
				}


				//std::cout << "Best elem for k = " << k << " is " << bestElemPerThread[thWithMaxElem] << std::endl;
				//std::cout << "best th is " << thWithMaxElem << std::endl;

				// 
				for (int th = 1; th < nprocs; ++th) {
					MPI_Send(&thWithMaxElem, 1, MPI_INT, th, 0, MPI_COMM_WORLD); // отправляем всем тредам максимальный элемент , а имеено индекс той нити где тот находжится
				}
				//std::cout << "Sended info " << std::endl;
				if (currentThread == thWithMaxElem) { // если это максимальный 
					int flag = 0;
					
					// 
					if (currentThread == 0) {
						for (int j = 0; j < n + 1; ++j) {
							double tmp = AB[k * (n + 1) + j];
							AB[k * (n + 1) + j] = AB[curBestRow * (n + 1) + j];
							AB[curBestRow * (n + 1) + j] = tmp;
						}
					}
				}
				else {

					// 
					if (thWithMaxElem == 0) {
						// send row with max
						MPI_Send(AB+curBestRow*(n+1), n+1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD);
						// recv row
						MPI_Recv(AB + curBestRow * (n + 1), n + 1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
						std::cout << "WRONG!!!!!!!!!!!!!!!!!!!" << std::endl;
						MPI_Finalize();
						return 1;
					}
					else {

						if (currentThread == 0) {// то есть мастер тред 1 проход но нет там максимальный элемент находится

							//std::cout << "0 <- " << thWithMaxElem << std::endl;
							MPI_Recv(R, n + 1, MPI_DOUBLE, thWithMaxElem, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
							//std::cout << "0 -> " << thWithMaxElem << std::endl;
							MPI_Send(AB + k * (n + 1), n + 1, MPI_DOUBLE, thWithMaxElem, 0, MPI_COMM_WORLD);
							// 
							for (int j = 0; j < n + 1; ++j) {
								AB[k * (n + 1) + j] = R[j];//вот тут не понятно, как будто мы должны не менять местами а просто з
							}// здесь вот не меняем разве%%
							//std::cout << "Done" << std::endl;
						}
						else {
							// nothing
						}
					}
				} // 


				// 
				if (rank == currentThread) {// 
					// if current, send then divide
					double Akk = AB[k * (n + 1)+k];
					//fprintf(pf, "Akk = %.4lf\n", Akk);
					for (int j = k; j < n + 1; ++j) {
						AB[k * (n + 1) + j] /= Akk;
					}

					for (int th = 0; th < nprocs; ++th) {
						if (th != rank) {
						//	fprintf(pf, "Sending row:\n");
							/*
							for (int j = 0; j < n + 1; ++j) {
							//	fprintf(pf, "%.5lf  ", AB[k * (n + 1) + j]);
							}
							fprintf(pf, "\n");
							fflush(pf);
							*/
							MPI_Send(AB + k * (n + 1), n + 1, MPI_DOUBLE, th, 0, MPI_COMM_WORLD);
						}
					}

					// local job

					for (int row = 0; row < rootrows; ++row) {
						if (row == k)continue;
						double Aik = AB[row * (n + 1)+k];
						for (int j = k; j < n + 1; ++j) {

							//fprintf(pf, "%.5lf -= %.5lf = ", AB[row * (n + 1) + j], AB[k * (n + 1) + j]);
							AB[row * (n + 1) + j] -= Aik * AB[k * (n + 1) + j];
							//fprintf(pf, "%.5lf\n", AB[row * (n + 1) + j]);
						}
					}
					


				}
				else {

					// not current
					// recv then subtract
					MPI_Recv(R, n + 1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (int row = 0; row < rootrows; ++row) {

						double Aik = AB[row * (n + 1)+k];
						for (int j = k; j < n + 1; ++j) {

							AB[row * (n + 1) + j] -= Aik * R[j];
						}
					}

				}

				
			}
			
			for (int i = 0; i < rootrows; ++i) {
				X[i] = AB[i * (n + 1) + n];
			}
			for (int th = 1; th < nprocs; ++th) { // счетчик
				// 
				MPI_Recv(R+th*nrows, nrows, MPI_DOUBLE, th, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // здесь он получает решения
				for (int i = 0; i < nrows; ++i) {
					X[rootrows+ (th-1) * nrows + i] = R[th * nrows + i];
				}
			}

			
			time_t tFinish = clock();

			for (int i = 0; i < n; ++i) {
				cout<<X[i]<<endl; 
				//fprintf(pf, "X[%d] = %.8lf\n", i, X[i]);
			}
			

			//fprintf(pf, "Time elapsed = %.3lf\n", 1.0 * (tFinish - tStart) / CLOCKS_PER_SEC);
			//fclose(pf);
			delete[] AB;
			delete[] X;
			delete[] B;
		}

	}// вот здесь заканчиваем с мастером
	else {
		// not master 
		n = stoi(argv[1]); 
		if (nprocs < n) {

			//char bufname[16] = {0}; 
			//snprintf(bufname, 16, "data#%d.txt", rank);
			
			int nrows = n / nprocs;
			int rootrows = n - (nprocs - 1) * nrows;// все то же самое что и для мастера делаем
			int rowStart = n - nrows * (nprocs - rank); // // начало в матице 

			

			printf("rank=%d  nrows=%d  rowStart=%d\n", rank, nrows, rowStart);

			double* AZ = new double[nrows * (n+1)]; // локальный массив(он состоит из тех строк которые мы передаем в родительский процесс)
			double* R = new double[n + 1];// заводим буффер 
			
			MPI_Recv(AZ, nrows * (n+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // тут наш процесс получает этот массив
			
			//FILE* pf = fopen(bufname, "w");
			/*
			for (int row = 0; row < nrows; ++row) {
				//fprintf(pf, "AB[%d] = ", rowStart + row);
				for (int col = 0; col < n+1; ++col) {
					//fprintf(pf, "%.5lf  ", AZ[row * (n+1) + col]);
				}
				fprintf(pf, "\n");
			}
			*/
			MPI_Recv(R, n + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE); // получаем 1 строку матрицы
			// first step;
			//fprintf(pf, "Received: \n");
			/*
			for (int j = 0; j < n + 1; ++j) {
				fprintf(pf, "%.5lf  ", R[j]);
			}
			fprintf(pf, "\n");
			*/
			
			for (int row = 0; row < nrows; ++row) { // вот здесб
				double Ai1 = AZ[row*(n+1)];
				for (int j = 0; j < n + 1; ++j) { // 
					AZ[row * (n + 1) + j] -= Ai1 * R[j];
				}
			} // выполняем 1 шаг 
			
			/*
			for (int row = 0; row < nrows; ++row) {
				fprintf(pf, "AB[%d] = ", rowStart + row);
				for (int col = 0; col < n + 1; ++col) {
					fprintf(pf, "%.5lf  ", AZ[row * (n + 1) + col]);
				}
				fprintf(pf, "\n");
			}
			*/
			for (int k = 1; k < n; ++k) {

				int currentThread = 0; // делаем то же самое что и в мастер потоке 
				if (k >= rootrows) {
					currentThread = 1 + (k - rootrows) / nrows;
				}
				//fprintf(pf, "\nk=%d  curThread=%d----------------------------------------\n", k, currentThread);
				/*
				for (int row = 0; row < nrows; ++row) {
					fprintf(pf, "AB[%d] = ", rowStart + row);
					for (int col = 0; col < n + 1; ++col) {
						fprintf(pf, "%.5lf  ", AZ[row * (n + 1) + col]);
					}
					fprintf(pf, "\n");
				}
				fflush(pf);
				*/
				
				
				// find max
				double bestAik = 0; //ищем максимальный 
				int bestRow = 0; // лучшая сторка
				
				if (rank != currentThread) {
					for (int i = 0; i < nrows; ++i) {
						if (fabs(AZ[i * (n + 1) + k]) > fabs(bestAik)) {
							bestAik = AZ[i * (n + 1) + k];
							bestRow = i;
						}
					}
				}
				else {
					int k2 = k - rootrows - (rank - 1) * nrows;
					bestAik = AZ[k2 * (n + 1) + k];
					bestRow = k2;
					for (int i = k2+1; i < nrows; ++i) {
						if (fabs(AZ[i * (n + 1) + k]) > fabs(bestAik)) {
							bestAik = AZ[i * (n + 1) + k];
							bestRow = i;
						}
					}
				}
				// 
				MPI_Send(&bestAik, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // отправляем мастер потоку 
				//std::cout << "rank " << rank << " sended " << bestAik << std::endl;


				// 
				int flag = 0;
				MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE); // принимаем максимальный эемент 
				

				//std::cout << "rank " << rank << "received flag " << flag << std::endl;

				if (flag == rank) { // если это наш элемент индекс, если быть точнее 

					// this thread has best row

					if (currentThread == rank) { // если одна из строк потока k-z

						// swap local
						
						int k2 = k - rootrows - (rank - 1) * nrows; // это локальный k 
					//	fprintf(pf, "rank %d swaps inside %d <-> %d\n", rank, bestRow, k2);
						if (k2 != bestRow) { //если локальный k не совпадает с лучшей строкой, то меняем местами 
							for (int j = 0; j < n + 1; ++j) {

								double tmp = AZ[k2 * (n + 1) + j];
								AZ[k2 * (n + 1) + j] = AZ[bestRow * (n + 1) + j];
								AZ[bestRow * (n + 1) + j] = tmp;
							} // локальный обммен если каррент ранг 
						}
						else { // иначе надобности их менять местами никакой нет 
						//	fprintf(pf, "no need to swap\n");
						}

					}
					else {
						
						MPI_Send(AZ + bestRow * (n + 1), n+1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD);
						MPI_Recv(AZ + bestRow * (n + 1), n + 1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // иначе отправляем, получаем 
					}
				}
				else {

					// not best, but maybe current

					if (currentThread == rank) {
						int k2 = k - rootrows - (rank - 1) * nrows;
						//fprintf(pf, "rank %d receive big row from %d\n", rank, flag);
						// recv then send
						MPI_Recv(R, n + 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // получаем от номера с максимальным элементом нашу максимальную строку 
						/*
						for (int j = 0; j < n + 1; ++j) {
							fprintf(pf, "%.5lf  ", R[j]);
						}
						fprintf(pf, "\n");
						fflush(pf);
						*/

						MPI_Send(AZ + k2 * (n + 1), n + 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD); // строку передаем строку текущую, нашу же мен]ем на максимальную 
						for (int j = 0; j < n + 1; ++j) {
							AZ[k2 * (n + 1) + j] = R[j]; // вот тут меняем нашу строку на максимальную 
						}
					}
					else {
						//not current, not best
					}
				} // на этом этапе мы закончили обмен строчками 
				
				if (rank == currentThread) {
					int k2 = k - rootrows - (rank - 1) * nrows;
		
					double Akk = AZ[k2 * (n + 1)+k];
					for (int j = 0; j < n + 1; ++j) {
						AZ[k2 * (n + 1) + j] /= Akk;
					} // далее главную строчку мы делим 
					// 
					for (int th = 0; th < nprocs; ++th) {
						if (th != rank) {
							//fprintf(pf, "%d -> %d\n", rank, th);
							MPI_Send(AZ + k2 * (n + 1), n + 1, MPI_DOUBLE, th, 0, MPI_COMM_WORLD);
						} // эту главную строчку мы отправляем нашим остальным нитям 
					}

					// local job
					for (int row = 0; row < nrows; ++row) { // тут в свою очередь 
						if (row == k2)continue;
						double Aik = AZ[row * (n + 1) + k];
						//fprintf(pf, "Aik = %.5lf\n", Aik);
						for (int j = 0; j < n + 1; ++j) {
							//fprintf(pf, "%.5lf -= %.5lf", AZ[row * (n + 1) + j], Aik * AZ[k2*(n+1)+j]);
							AZ[row * (n + 1) + j] -= Aik * AZ[k2 * (n + 1) + j]; // вот тут наверно не сильно оптимально получается потому что мы будем вычитать нули, но короче в рамках своего куска проделываем нужную работу 
							//fprintf(pf, " = %.5lf\n", AZ[row * (n + 1) + j]);
						}
						//fprintf(pf, "\n");
					}
					//fflush(pf);

				}
				else {

					// not current
					// recv then subtract
					MPI_Recv(R, n + 1, MPI_DOUBLE, currentThread, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // иначе получаем эту максимальную строку 
					//fprintf(pf, "Received row:\n");
					/*
					for (int j = 0; j < n + 1; ++j) {
						fprintf(pf, "%.4lf  ", R[j]);
					}
					*/
					//fprintf(pf, "\n");
					//fflush(pf);
					for (int row = 0; row < nrows; ++row) {
						
						double Aik = AZ[row * (n + 1)+k];
						for (int j = 0; j < n + 1; ++j) {
							//fprintf(pf, "%.4lf -= %.4lf", AZ[row * (n + 1) + j], Aik * R[j]);
							AZ[row * (n + 1) + j] -= Aik * R[j]; // занимаемся тем же самвм
							//fprintf(pf, " = %.4lf\n", AZ[row * (n + 1) + j]);
						}
						//fprintf(pf, "\n");
					}
					//fflush(pf);


				}

				
			}

			for (int i = 0; i < nrows; ++i) {
				R[i] = AZ[i * (n + 1) + n];
			} //здесь записываем ответ конкретно для нашей части и отправляем в масер 
			
			MPI_Send(R, nrows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

			
			delete[] AZ;
			delete[] R;
		}

	}



	MPI_Finalize();
	return 0;
}