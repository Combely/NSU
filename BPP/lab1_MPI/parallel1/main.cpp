#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <mpi.h>

struct Matrix {
	double *matrix = NULL;
	int rows = 0;
	int columns = 0;
};

int N, rank, num_of_procs, chunk_start, chunk_end;
int *mpi_sizes, *mpi_displs;
Matrix buffer, buffer_complete;

Matrix MatrixInit(int rows, int columns);
void MatrixCheckSizes(Matrix& input, int exp_rows, int exp_columns);
void MatrixFree(Matrix& input);
void MatrixFillA(Matrix& A);
void MatrixFillRand(Matrix& input);
void MatrixFillZero(Matrix& zero);
void MatrixToMatrixCopy(Matrix& source, Matrix& dest);
void MatrixToMatrixAdd(Matrix& left, Matrix& right, Matrix& res);
void MatrixFromMatrixSub(Matrix& left, Matrix& right, Matrix& res);
void MatrixOnMatrixMult(Matrix& left, Matrix& right, Matrix& res);
void MatrixOnScalarMult(Matrix& input, double scalar, Matrix& res);
double VectorScalarMult(Matrix& lvect, Matrix& rvect);
double GetSquaredNorm(Matrix& input);
Matrix GetR0(Matrix& b, Matrix& A, Matrix& x);
Matrix GetZ0(Matrix& r);
double UpdateAlpha(Matrix& r, Matrix& A, Matrix& z);
void UpdateX(Matrix& x, double alpha, Matrix& z);
void UpdateR(Matrix& r, double alpha, Matrix& A, Matrix& z);
double UpdateBeta(Matrix& r, double prev_scalar_r);
void UpdateZ(Matrix& r, double beta, Matrix& z);
void PrepSizesAndDispls(int** sizes, int** displs, int num_of_procs, int raws, int columns);
void PrepLocalIndexes();

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	double start, end;
	int local_A_raws;

	if (!rank) {
		start = MPI_Wtime();
		if (argc < 2) {
			N = 1000;
		}
		else {
			N = atoi(argv[1]);
		}
		Matrix A = MatrixInit(N, N);
		MatrixFillA(A);
	}
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	buffer = MatrixInit(N, 1);
	buffer_complete = MatrixInit(N, 1);
	Matrix x = MatrixInit(N, 1);
	MatrixFillZero(x);
	Matrix b = MatrixInit(N, 1);
	MatrixFillRand(b);

	PrepSizesAndDispls(&mpi_sizes, &mpi_displs, num_of_procs, N, N);
	local_A_raws = mpi_sizes[rank] / N;
	Matrix local_A = MatrixInit(local_A_raws, N);
	MPI_Scatterv(A.matrix, mpi_sizes, mpi_displs, MPI_DOUBLE,
        local_A.matrix, local_A_raws, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	PrepLocalIndexes();

	Matrix r = GetR0(b, local_A, x);
	Matrix z = GetZ0(r);
	
	int res_criteria_cter = 0, cycles_cter = 0;
	double squared_norm_b = GetSquaredNorm(b);
	double squared_norm_r, alpha, beta;
	double epsilon = (1e-5) * (1e-5) * squared_norm_b;
	
	while (res_criteria_cter != 5) {
		cycles_cter++;
		if (cycles_cter > 10000) {
				fprintf(stderr, "Unable to get result\n");
				exit(EXIT_FAILURE);
		}
		if ((squared_norm_r = GetSquaredNorm(r)) < epsilon) {
			res_criteria_cter++;
		}
		else {
			res_criteria_cter = 0;
		}
		alpha = UpdateAlpha(r, local_A, z);
		UpdateX(x, alpha, z);
		UpdateR(r, alpha, local_A, z);
		beta = UpdateBeta(r, squared_norm_r);
		UpdateZ(r, beta, z);
	}
	
	if (!rank) {
		end = MPI_Wtime();
		printf("Version: Parallel v.1\n Matrix size: %d\n Consumed time: %lf\n", N, end - start);
		MatrixFree(A);
	}
	MatrixFree(local_A);
	MatrixFree(x);
	MatrixFree(b);
	MatrixFree(r);
	MatrixFree(z);
	MatrixFree(buffer);
	MatrixFree(buffer_complete);
	MPI_Finalize();
	return 0;
}

Matrix MatrixInit(int rows, int columns) {
	Matrix unfilled;
	unfilled.rows = rows;
	unfilled.columns = columns;
	unfilled.matrix = (double*)malloc(rows * columns * sizeof(double));
	if (unfilled.matrix == NULL) {
		printf("Memory is not allocated.\n");
		exit(0);
    }
	return unfilled;
}

void MatrixCheckSizes(Matrix& input, int exp_rows, int exp_columns) {
	if (input.rows != exp_rows || input.columns != exp_columns) {
		MatrixFree(input);
		input = MatrixInit(exp_rows, exp_columns);
	}
}

void MatrixFree(Matrix& input) {
	if (&input) {
		free(input.matrix);
	}
}

void MatrixFillA(Matrix& A) {
	srand(100); 
	double weight = 100;
	double* A_matrix = A.matrix;
	int rows = A.rows;
	int columns = A.columns;
	for (int row = 0; row < rows; row++) {
		for (int column = 0; column <= row; column++) {
			A_matrix[row * columns + column] = rand() % 11;
      		if (row == column) {
				A_matrix[row * columns + column] += weight;
				continue;
			}
      		A_matrix[column * columns + row] = A_matrix[row * columns + column];
		} 
	}
}

void MatrixFillRand(Matrix& input) {
	srand(59);
	int size = input.rows * input.columns;
	for (int i = 0; i < size; i++) {
		input.matrix[i] = rand() % 11;
	}
}

void MatrixFillZero(Matrix& zero) {	
	memset(zero.matrix, 0, sizeof(double) * zero.columns * zero.rows);
}

void MatrixToMatrixCopy(Matrix& source, Matrix& dest) {
	int size = source.rows * source.columns;
	memcpy(dest.matrix, source.matrix, sizeof(double) * size);
}

void MatrixToMatrixAdd(Matrix& left, Matrix& right, Matrix& res) {
	int size = left.rows * left.columns;
	double* left_matrix = left.matrix;
	double* right_matrix = right.matrix;
	double* res_matrix = res.matrix;
	MatrixCheckSizes(res, right.rows, right.columns);
	for (int i = chunk_start; i <= chunk_end; i++) {
		res_matrix[i] = left_matrix[i] + right_matrix[i];
	}
	MPI_Allgatherv(res_matrix + chunk_start, mpi_sizes[rank],
        MPI_DOUBLE, res_matrix, mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

void MatrixFromMatrixSub(Matrix& left, Matrix& right, Matrix& res) {
	int size = left.rows * left.columns;
	double* left_matrix = left.matrix;
	double* right_matrix = right.matrix;
	double* res_matrix = res.matrix;
	MatrixCheckSizes(res, right.rows, right.columns);
	for (int i = chunk_start; i <= chunk_end; i++) {
		res_matrix[i] = left_matrix[i] - right_matrix[i];
	}
	MPI_Allgatherv(res_matrix + chunk_start, mpi_sizes[rank],
        MPI_DOUBLE, res_matrix, mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

void MatrixOnMatrixMult(Matrix& left, Matrix& right, Matrix& res) {
	int r1 = left.rows, c1 = left.columns; 
	int r2 = right.rows, c2 = right.columns;
	double* left_matrix = left.matrix;
	double* right_matrix = right.matrix;
	double* res_matrix = res.matrix;
	MatrixCheckSizes(res, r1, c2);
	for (int i = 0; i < r1; i++) {
		for (int j = 0; j < c2; j++) {
			res_matrix[i * c2 + j] = 0;
			for (int k = 0; k < r2; k++) {
				res_matrix[i * c2 + j] += left_matrix[i * c1 + k] * right_matrix[k * c2 + j];
			}
		}
	}

}

void MatrixOnScalarMult(Matrix& input, double scalar, Matrix& res) {
	int size = input.rows * input.columns;
	double* input_matrix = input.matrix;
	double* res_matrix = res.matrix;
	MatrixCheckSizes(res, input.rows, input.columns);
	for (int i = chunk_start; i <= chunk_end; i++) {
		res_matrix[i] = input_matrix[i] * scalar;
	}
	MPI_Allgatherv(res_matrix + chunk_start, mpi_sizes[rank],
        MPI_DOUBLE, res_matrix, mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
}


double VectorScalarMult(Matrix& lvect, Matrix& rvect) {
	double res = 0, part_sum = 0;
	double* lvect_data = lvect.matrix;
	double* rvect_data = rvect.matrix;
	int size = lvect.rows;
	for (int i = chunk_start; i < chunk_end; i++) {
		part_sum += lvect_data[i] * rvect_data[i];
	}
	MPI_Allreduce(&part_sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return res;
}

double GetSquaredNorm(Matrix& input) {
	return VectorScalarMult(input, input);
}

Matrix GetR0(Matrix& b, Matrix& A_chunk, Matrix& x) {
	Matrix r0 = MatrixInit(N, 1);
	MatrixOnMatrixMult(A_chunk, x, buffer);
	MPI_Allgatherv(buffer.matrix, mpi_sizes[rank], MPI_DOUBLE, r0.matrix, 
		mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
	MatrixFromMatrixSub(b, r0, r0);
	return r0;
}

Matrix GetZ0(Matrix& r) {
	Matrix z = MatrixInit(N, 1);
	MatrixToMatrixCopy(r, z);
	return z;
}

double UpdateAlpha(Matrix& r, Matrix& A_chunk, Matrix& z) {
	MatrixOnMatrixMult(A_chunk, z, buffer);
	MPI_Allgatherv(buffer.matrix, mpi_sizes[rank], MPI_DOUBLE,
        buffer_complete.matrix, mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
	double numerator = VectorScalarMult(r, r);
	double denominator = VectorScalarMult(buffer_complete, z);
	return (numerator / denominator);
}

void UpdateX(Matrix& x, double alpha, Matrix& z) {
	MatrixOnScalarMult(z, alpha, buffer);
	MatrixToMatrixAdd(x, buffer, x);
}
void UpdateR(Matrix& r, double alpha, Matrix& A_chunk, Matrix& z) {
	MatrixOnMatrixMult(A_chunk, z, buffer);
	MPI_Allgatherv(buffer.matrix, mpi_sizes[rank], MPI_DOUBLE,
        buffer_complete.matrix, mpi_sizes, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
	MatrixOnScalarMult(buffer_complete, alpha, buffer_complete);
	MatrixFromMatrixSub(r, buffer_complete, r);
}

double UpdateBeta(Matrix& r, double prev_scalar_r) {
	double curr_scalar_r = VectorScalarMult(r, r);
	return (curr_scalar_r / prev_scalar_r);
}

void UpdateZ(Matrix& r, double beta, Matrix& z) {
	MatrixOnScalarMult(z, beta, z);
	MatrixToMatrixAdd(r, z, z);
}

void PrepSizesAndDispls(int** sizes, int** displs, int num_of_procs, int raws, int columns) {
	int main_chunk = rows / num_of_procs;
  	int remainder = rows % num_of_procs;
	int offset = 0;
	int* pter_sizes = *sizes, *pter_displs = *displs;
  	pter_sizes = (int*)malloc(sizeof(int) * num_of_procs);
  	pter_displs = (int*)malloc(sizeof(int) * num_of_procs);
  	for (int i = 0; i < num_of_proc; i++) {
		pter_sizes[i] = main_chunk * columns;
    	pter_displs[i] = offset;
    	offset += main_chunk * columns;
    	if (remainder) { //each process takes no more than 1
      		pter_sizes[i] += columns;
      		offset += columns;
     		remainder--;
    	}
  	}
}

void PrepLocalIndexes() {
	for (int i = 0; i < num_of_procs; i++) {
		mpi_displs[i] /= N;
		mpi_sizes[i] /= N;
	}
	chunk_start = mpi_displs[rank];
	chunk_end = mpi_displs[rank] + mpi_sizes[rank] - 1;
}
