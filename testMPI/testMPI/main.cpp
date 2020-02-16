#include <stdio.h>
#include "mpi.h"
#include <sys\timeb.h>
#include <SDKDDKVer.h>
#include <cstdlib>
#include <tchar.h>
#include <time.h>
#define _CRT_SECURE_NO_WARNINGS

int main(int argc, char ** argv)
{
	srand(time(0));

	int msize = 10000; // всего строк в матрице
	int rl = 0; // row length (вводим с клавиатуры)
	int** a = NULL; // matrix
	int* v = NULL; // vector
	int myrank, ranksize;
	MPI_Status status;
	double timer = 0;
	int rcnt = 0; //сколько строк обрабатывает каждый процесс

	MPI_Init(&argc, &argv);		/* initialize MPI system */
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);		/* my place in MPI system */
	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);	/* size of MPI system */

	if (myrank == 0)	/* I am the master */
	{
		printf("Starting lab1...\n");
		scanf_s("%d", &rl);
		printf("\n");

		a = new int*[msize];
		v = new int[rl];
		for (int i = 0; i < msize; i++) {
			a[i] = new int[rl];
			for (int j = 0; j < rl; j++) {
				a[i][j] = i + 1;
			}
		}
		for (int j = 0; j < 5; j++) {
			v[j] = 1;
		}
		for (int j = 5; j < rl; j++) {
			v[j] = 0;
		}

		timer -= MPI_Wtime();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&rl, 1, MPI_INT, 0, MPI_COMM_WORLD); //послали всем длину строки/вектора
	MPI_Barrier(MPI_COMM_WORLD);

	rcnt = msize / ranksize; //сколько строк обрабатывает каждый процесс
	if (myrank != 0) {
		v = new int[rl];
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(v, rl, MPI_INT, 0, MPI_COMM_WORLD); //послали сам вектор
	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0) {
		int i = rcnt + msize % ranksize; //root забирает себе один кусок + все лишние строки
		for (int k = 1; k < ranksize; k++) { //идем по процессам
			for (int j = 0; j < rcnt; j++) {
				MPI_Send(a[i], rl, MPI_INT, k, 100 + j, MPI_COMM_WORLD);//отправили строку
				i++;
			}
		}
	}
	else {
		a = new int*[rcnt];
		for (int i = 0; i < rcnt; i++) {
			a[i] = new int[rl];
			MPI_Recv(a[i], rl, MPI_INT, 0, 100 + i, MPI_COMM_WORLD, &status);//если мы не root, то приняли строку
		}
	}

	if (myrank == 0) printf("\n Had sent data - Execution time: %.3lf sec\n", timer + MPI_Wtime()); //сколько времени ушло на посылку

	long* ans = NULL;
	if (myrank == 0)	// I am the master, collect results, add up, and print results
	{
		ans = new long[msize];
		for (int i = 0; i < rcnt + msize % ranksize; i++) {
			ans[i] = 0;
			for (int j = 0; j < rl; j++) {
				ans[i] += a[i][j] * v[j];
			}
		} //root посчитал свою часть

		for (int i = 1; i < ranksize; i++)//собираем данные от всех процессов
		{
			long* buf = new long[rcnt];
			MPI_Recv(buf, rcnt, MPI_LONG, i, i, MPI_COMM_WORLD, &status);

			for (int j = 0; j < rcnt; j++) {
				ans[j + msize % ranksize + rcnt * i] = buf[j];
			}
		}

		timer += MPI_Wtime();
		for (int i = 0; i < 10; i++) {
			printf("%ld ", ans[i]);
		}
		printf("\n%d tasks used - Execution time: %.3lf sec\n", ranksize, timer);
	}
	else {	// I am a slave
			// send my result back to master 
		ans = new long[rcnt];
		for (int i = 0; i < rcnt; i++) {
			ans[i] = 0;
			for (int j = 0; j < rl; j++) {
				ans[i] += a[i][j] * v[j];
			}
		}
		MPI_Send(ans, rcnt, MPI_LONG, 0, myrank, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}

