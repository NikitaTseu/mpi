#define _CRT_SECURE_NO_WARNINGS

#include<mpi.h>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char ** argv)
{
	srand(17);

	MPI_Status status;
	int myrank, ranksize;
	double t1, t2;
	double **a = NULL;
	int* imap = NULL;
	int n = 10;
	std::string filename = "input_00.txt";
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	/*if (myrank == 0)
	{
		std::ofstream fout("input_00.txt");
		int size = 1000;

		fout << size << std::endl;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++) {
				fout << rand() % 10000 << " ";
			}
			fout << std::endl;
		}
		fout.close();
	}*/

	if (myrank == 0) {
		std::cin >> filename;
		std::ifstream fin(filename);
		fin >> n;
		fin.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	imap = new int[n];
	for (int i = 0; i < n; i++) {
		imap[i] = i % ranksize;
	}

	a = new double*[n];
	for (int i = 0; i < n; i++) {
		a[i] = new double[n];
	}

	if (myrank == 0)
	{
		std::ifstream fin(filename);
		fin >> n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				fin >> a[i][j];
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < n; i++) {
		MPI_Bcast(a[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	t1 = MPI_Wtime();

	for (int k = 0; k < n - 1; k++)
	{
		for (int i = k + 1; i < n; i++) {
			if (imap[i] == myrank) {
				a[i][k] /= a[k][k];
				for (int j = k + 1; j < n; j++) {
					a[i][j] -= a[k][j] * a[i][k];
				}

			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(a[k+1], n, MPI_DOUBLE, imap[k+1], MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();

	if (myrank == 0)
	{
		std::ofstream fout("output_10.txt");

		fout << "L:" << std::endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++) {
				fout << a[i][j] << " ";
			}	
			fout << "1 ";
			for (int j = i + 1; j < n; j++) {
				fout << "0 ";
			}
			fout << std::endl;
		}

		fout << std::endl << "U:" << std::endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++) {
				fout << "0 ";
			}
			for (int j = i; j < n; j++) {
				fout << a[i][j] << " ";
			}
			fout << std::endl;
		}
	}

	if (myrank == 0) {
		printf("\n%d tasks used - Execution time: %.3lf sec\n", ranksize, t2 - t1);
	}

	MPI_Finalize();
	return 0;
}

