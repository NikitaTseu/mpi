#define _CRT_SECURE_NO_WARNINGS

#include<mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

int N = 101;

double h = 1.0 / N;

double f(double x, double y) {return x * y;}

double f1(double y) { return sin(y); }
double f2(double y) { return y; }
double f3(double x) { return x * x * x; }
double f4(double x) { return x * x; }

double value(double left, double right, double up, double down, double f) {
	return 0.25 * (left + right + up + down - h * h * f);
}

int main(int argc, char** argv) {

	MPI_Status status;
	int myrank, ranksize;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int n = ((N + 1) - 2) / ranksize;

	double* grid = new double[N + 1];
	for (int i = 0; i < N + 1; i++) {
		grid[i] = i * h;
	}

	double** a = new double*[n];
	double** b = new double*[n];

	for (int i = 0; i < n; i++) {
		a[i] = new double[N + 1];
		b[i] = new double[N + 1];

		for (int j = 0; j < N + 1; j++) {
			a[i][j] = 0.0;
			b[i][j] = 0.0;
		}
	}

	double* up = new double[N + 1];
	double* down = new double[N + 1];
	for (int j = 0; j < N + 1; j++) {
		up[j] = 0.0;
		down[j] = 0.0;
	}
	
	if (myrank == 0) {
		for (int j = 0; j < N + 1; j++) {
			up[j] = f4(grid[j]);
		}
	}

	if (myrank == ranksize - 1) {
		for (int j = 0; j < N + 1; j++) {
			down[j] = f3(grid[j]);
		}
	}

	for (int i = 0; i < n; i++) {
		a[i][0] = f1(grid[myrank*n + i + 1]);
		a[i][N] = f2(grid[myrank*n + i + 1]);
	}

	// iterations
	double time = MPI_Wtime();

	for (int k = 0; k < 1000; k++) {
		for (int j = 1; j < N; j++) {
			b[0][j] = value(a[0][j - 1], a[0][j + 1], up[j], a[1][j], f(grid[j], grid[myrank*n + 1]));
		}

		for (int j = 1; j < N; j++) {
			b[n-1][j] = value(a[n - 1][j - 1], a[n - 1][j + 1], a[n - 2][j], down[j], f(grid[j], grid[myrank*n + (n-1) + 1]));
		}

		for (int i = 1; i < n-1; i++) {
			for (int j = 1; j < N; j++) {
				b[i][j] = value(a[i][j - 1], a[i][j + 1], a[i - 1][j], a[i + 1][j], f(grid[j], grid[myrank*n + i + 1]));
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 1; j < N; j++) {
				a[i][j] = b[i][j];
 			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (ranksize > 1) {
			if (myrank == 0) {
				MPI_Send(a[n - 1], N + 1, MPI_DOUBLE, 1, 3000, MPI_COMM_WORLD);
				MPI_Recv(down, N + 1, MPI_DOUBLE, 1, 1001, MPI_COMM_WORLD, &status);
			}

			if (myrank == ranksize - 1) {
				MPI_Send(a[0], N + 1, MPI_DOUBLE, ranksize - 2, 1000 + ranksize - 1, MPI_COMM_WORLD);
				MPI_Recv(up, N + 1, MPI_DOUBLE, ranksize - 2, 3000 + ranksize - 2, MPI_COMM_WORLD, &status);
			}

			if ((myrank != 0) && (myrank != ranksize - 1)) {
				MPI_Send(a[0], N + 1, MPI_DOUBLE, myrank - 1, 1000 + myrank, MPI_COMM_WORLD);
				MPI_Send(a[n - 1], N + 1, MPI_DOUBLE, myrank + 1, 3000 + myrank, MPI_COMM_WORLD);
				MPI_Recv(up, N + 1, MPI_DOUBLE, myrank - 1, 3000 + myrank - 1, MPI_COMM_WORLD, &status);
				MPI_Recv(down, N + 1, MPI_DOUBLE, myrank + 1, 1000 + myrank + 1, MPI_COMM_WORLD, &status);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0) {
		std::cout << "Time: " << MPI_Wtime() - time << std::endl;
	}

	std::ofstream out("matrix.txt");

	double** ans = new double*[N + 1];

	if (myrank == 0) {
		for (int i = 0; i < N + 1; i++) {
			ans[i] = new double[N + 1];
		}

		for (int j = 0; j < N + 1; j++) {
			ans[0][j] = up[j];
			ans[N][j] = f3(grid[j]);
		}

		for (int i = 1; i < n + 1; i++) {
			for (int j = 0; j < N + 1; j++) {
				ans[i][j] = a[i - 1][j];
			}
		}

		for (int i = n + 1; i < N; i++) {
			MPI_Recv(ans[i], N + 1, MPI_DOUBLE, (i-1)/n, i, MPI_COMM_WORLD, &status);
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			MPI_Send(a[i], N + 1, MPI_DOUBLE, 0, myrank*n + 1 + i, MPI_COMM_WORLD);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0) {
		for (int i = 0; i < N + 1; i++) {
			for (int j = 0; j < N + 1; j++)
				out << ans[i][j] << " ";
			out << std::endl;
		}
	}

	MPI_Finalize();
	return 0;
}