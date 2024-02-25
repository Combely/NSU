#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <omp.h>
#include <thread>
struct Matrix {
	double *matrix = NULL;
	int rows = 0;
	int columns = 0;
};

int N;
double result;
Matrix buffer, shared;

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
void VectorScalarMult(Matrix& lvect, Matrix& rvect, volatile double& res);
void GetSquaredNorm(Matrix& input, volatile double& squared_res);
Matrix GetR0(Matrix& b, Matrix& A, Matrix& x);
Matrix GetZ0(Matrix& r);
double UpdateAlpha(Matrix& r, Matrix& A, Matrix& z);
void UpdateX(Matrix& x, double alpha, Matrix& z);
void UpdateR(Matrix& r, double alpha, Matrix& A, Matrix& z);
double UpdateBeta(Matrix& r, double prev_scalar_r);
void UpdateZ(Matrix& r, double beta, Matrix& z);

int main(int argc, char** argv) {
	if (argc < 3) {
		omp_set_num_threads(2);
		if (argc == 2) {
			N = atoi(argv[1]);
		}
		else {
			N = 1000;
		}
	}
	else {
		N = atoi(argv[1]);
    		omp_set_num_threads(atoi(argv[2]));
	}
	
	double start = omp_get_wtime();
	Matrix A, x, b, r, z;
	double vec_res = 0;

	#pragma omp parallel
	{

		#pragma omp single 
		{
			buffer = MatrixInit(N, 1);
			A = MatrixInit(N, N);
			MatrixFillA(A);
			x = MatrixInit(N, 1);
			MatrixFillZero(x);
			b = MatrixInit(N, 1);
			MatrixFillRand(b);
		}
		#pragma barrier
		Matrix local = GetR0(b, A, x);
		#pragma omp single nowait
		{
			r.rows = local.rows;
			r.columns = local.columns;
			r.matrix = local.matrix;
			shared.rows = 0;
			shared.columns = 0;
			shared.matrix = NULL;
			z = GetZ0(r);
		}
		int res_criteria_cter = 0, cycles_cter = 0;
		VectorScalarMult(b, b, vec_res);
		#pragma omp barrier
		double squared_norm_b = result;
		double squared_norm_r, alpha, beta;
		#pragma omp barrier
		double epsilon = (1e-5) * (1e-5) * squared_norm_b;
		printf("Thread %d, squared norm b = %lf, res = %lf\n", omp_get_thread_num(), squared_norm_b, result);
		#pragma omp barrier
		while (res_criteria_cter != 5) {
			cycles_cter++;
			if (cycles_cter > 10000) {
					fprintf(stderr, "Unable to get result\n");
					exit(EXIT_FAILURE);
			}
			GetSquaredNorm(r, squared_norm_r);
			if (squared_norm_r < epsilon) {
				res_criteria_cter++;
			}
			else {
				res_criteria_cter = 0;
			}

			alpha = UpdateAlpha(r, A, z);
			UpdateX(x, alpha, z);
			UpdateR(r, alpha, A, z);
			beta = UpdateBeta(r, squared_norm_r);
			UpdateZ(r, beta, z);
		}
	}
	
	
	double end = omp_get_wtime();
	
	printf("Version: Basic\n Matrix size: %u\n Consumed time: %lf\n", N, end - start);
	MatrixFree(A);
	MatrixFree(x);
	MatrixFree(b);
	MatrixFree(r);
	MatrixFree(z);
	MatrixFree(buffer);
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
	#pragma omp single
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
	#pragma omp for 
	for (int i = 0; i < size; i++) {
		res_matrix[i] = left_matrix[i] + right_matrix[i];
	}
}

void MatrixFromMatrixSub(Matrix& left, Matrix& right, Matrix& res) {
	int size = left.rows * left.columns;
	double* left_matrix = left.matrix;
	double* right_matrix = right.matrix;
	double* res_matrix = res.matrix;
	#pragma omp single
	MatrixCheckSizes(res, right.rows, right.columns);
	#pragma omp for 
	for (int i = 0; i < size; i++) {
		res_matrix[i] = left_matrix[i] - right_matrix[i];
	}
}

void MatrixOnMatrixMult(Matrix& left, Matrix& right, Matrix& res) {
	int r1 = left.rows, c1 = left.columns; 
	int r2 = right.rows, c2 = right.columns;
	double* left_matrix = left.matrix;
	double* right_matrix = right.matrix;
	double* res_matrix = res.matrix;
	#pragma omp single
	MatrixCheckSizes(res, r1, c2);
	#pragma omp for
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
	#pragma omp single
	MatrixCheckSizes(res, input.rows, input.columns);
	#pragma omp for
	for (int i = 0; i < size; i++) {
		res_matrix[i] = input_matrix[i] * scalar;
	}
}


void VectorScalarMult(Matrix& lvect, Matrix& rvect, volatile double& vec_res) {
	#pragma omp barrier
	#pragma omp single
	result = 0;
	double* lvect_data = lvect.matrix;
	double* rvect_data = rvect.matrix;
	int size = lvect.rows;
	#pragma omp barrier
	#pragma omp for reduction(+ : result)
	for (int i = 0; i < size; i++) {
		result += lvect_data[i] * rvect_data[i];
	}
	//#pragma omp barrier
}

void GetSquaredNorm(Matrix& input, volatile double &squared_norm) {
	#pragma omp barrier
	VectorScalarMult(input, input, squared_norm);
}

Matrix GetR0(Matrix& b, Matrix& A, Matrix& x) {
	#pragma omp single
	shared = MatrixInit(N, 1);
	MatrixOnMatrixMult(A, x, buffer);
	MatrixFromMatrixSub(b, buffer, shared);
	return shared;
}

Matrix GetZ0(Matrix& r) {
	Matrix z = MatrixInit(N, 1);
	MatrixToMatrixCopy(r, z);
	return z;
}

double UpdateAlpha(Matrix& r, Matrix& A, Matrix& z) {
	MatrixOnMatrixMult(A, z, buffer);
	double numerator, denominator;
	VectorScalarMult(r, r, numerator);
	VectorScalarMult(buffer, z, denominator);
	return (numerator / denominator);
}

void UpdateX(Matrix& x, double alpha, Matrix& z) {
	MatrixOnScalarMult(z, alpha, buffer);
	MatrixToMatrixAdd(x, buffer, x);
}
void UpdateR(Matrix& r, double alpha, Matrix& A, Matrix& z) {
	MatrixOnMatrixMult(A, z, buffer);
	MatrixOnScalarMult(buffer, alpha, buffer);
	MatrixFromMatrixSub(r, buffer, r);
}

double UpdateBeta(Matrix& r, double prev_scalar_r) {
	double curr_scalar_r;
	VectorScalarMult(r, r, curr_scalar_r);
	return (curr_scalar_r / prev_scalar_r);
}

void UpdateZ(Matrix& r, double beta, Matrix& z) {
	MatrixOnScalarMult(z, beta, z);
	MatrixToMatrixAdd(r, z, z);
}
