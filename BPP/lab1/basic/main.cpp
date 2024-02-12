#include <memory.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>

struct Matrix {
	double *matrix = NULL;
	int rows = 0;
	int columns = 0;
};

int N;

Matrix MatrixInit(int rows, int columns) {
	Matrix unfilled;
	unfilled.rows = rows;
	unfilled.columns = columns;
	unfilled.matrix = new double[rows * columns];
	return unfilled;
}

void MatrixFree(Matrix& input) {
	if (&input) {
		delete[] input.matrix;
	}
}

void MatrixFillA(Matrix& A) {
	double weight = 100;
	int rows = A.rows;
	int columns = A.columns;
	for (int row = 0; row < rows; row++) {
		for (int column = 0; column < columns; column++) {
			A.matrix[row * rows + column] = rand() % 11;
			if (row = column) {
				A.matrix[row * rows + column] += weight;
			}
		}
	}
}

void MatrixFillRand(Matrix& input) {
	int size = input.rows * input.columns;
	for (int i = 0; i < size; i++) {
		input.matrix[i] = rand() % 11;
	}
}

void MatrixFillZero(Matrix& zero) {	
	int size = zero.rows * zero.columns;
	for (int i = 0; i < size; i++) {
		zero.matrix[i] = 0;
	}
}

void MatrixOnMatrixMult(Matrix& left, Matrix& right, Matrix& res) {
	int r1 = left.rows, c1 = left.columns;
	int r2 = right.rows, c2 = right.columns;
	double* leftmatrix = left.matrix;
	double* rightmatrix = right.matrix;
	double* resmatrix = res.matrix;
	for (int i = 0; i < r1; i++) {
		for (int j = 0; j < c2; j++) {
			res.matrix[i * r1 + j] = 0;
			for (int k = 0; k < r2; k++) {
				resmatrix[i * r1 + j] += leftmatrix[i * r1 + k] * rightmatrix[k * r2 + j];
			}
		}
	}

}

void MatrixFromMatrixSub(Matrix& left, Matrix& right, Matrix& res) {
	int size = left.rows + left.columns;
	double* leftmatrix = left.matrix;
	double* rightmatrix = right.matrix;
	double* resmatrix = res.matrix;
	for (int i = 0; i < size; i++) {
		resmatrix[i] = leftmatrix[i] - rightmatrix[i];
	}
}

void MatrixOnScalarMult(Matrix& input, double scalar, Matrix& res) {
	int size = input.rows + input.columns;
	double* input_matrix = input.matrix;
	double* res_matrix = res.matrix;
	for (int i = 0; i < size; i++) {
		res_matrix[i] = input_matrix[i] * scalar;
	}
}


double VectorScalarMult(Matrix& lvect, Matrix& rvect) {
	double res = 0;
	double* lvectdata = lvect.matrix;
	double* rvectdata = rvect.matrix;
	int size = lvect.rows;
	for (int i = 0; i < size; i++) {
		res += lvectdata[i] * rvectdata[i];
	}
	return res;
}

void UpdateY(Matrix& y, Matrix& A, Matrix& x, Matrix& b) {
	MatrixOnMatrixMult(A, x, y);
	MatrixFromMatrixSub(y, b, y);
}

double GetTau(Matrix& A, Matrix& y) {
	double tau, numerator, denominator;
	Matrix buffer = MatrixInit(N, 1);
	MatrixOnMatrixMult(A, y, buffer);
	numerator = VectorScalarMult(y, buffer);
	denominator = VectorScalarMult(buffer, buffer);
	tau = numerator / denominator;
	MatrixFree(buffer);
	return tau;
}

void UpdateX(Matrix& x, Matrix& y, double tau) {
	Matrix buffer = MatrixInit(N, 1);
	MatrixOnScalarMult(y, tau, buffer);
	MatrixFromMatrixSub(x, buffer, x);
	MatrixFree(buffer);
}

double GetSquaredNorm(Matrix& input) {
	double res = 0;
	double* vector = input.matrix;
	int size = input.rows;
	for (int i = 0; i < size; i++) {
		res += vector[i] * vector[i];
	}
	return res;
}

int main(int argc, char** argv) {
	if (argc < 2) {
		N = 1000;
	}
	else {
		N = atoi(argv[1]);
	}
	time_t start = time(NULL);
	Matrix A = MatrixInit(N, N);
	MatrixFillA(A);
	Matrix x = MatrixInit(N, 1);
	MatrixFillZero(x);
	Matrix b = MatrixInit(N, 1);
	MatrixFillRand(b);
	Matrix y = MatrixInit(N, 1);
	UpdateY(y, A, x, b);
	int res_criteria_cter = 0, cycles_cter = 0;
	double squared_norm_b = GetSquaredNorm(b);
	double squared_norm_y, tau;
	double epsilon = (1e-5) * (1e-5) * squared_norm_b;
	while (res_criteria_cter != 5) {
		cycles_cter++;
                if (cycles_cter > 10000) {
                        fprintf(stderr, "Unable to get result");
                        exit(EXIT_FAILURE);
                }
		if (squared_norm_y = GetSquaredNorm(y) < epsilon) {
			res_criteria_cter++;
		}
		else {
			res_criteria_cter = 0;
		}

		tau = GetTau(A, y);
		UpdateX(x, y, tau);
		UpdateY(y, A, x, b);
	}
	time_t end = time(NULL);
	printf("Version: Basic\n Matrix size: %u\n Consumed time: %ld", N, end - start);
	MatrixFree(A);
	MatrixFree(x);
	MatrixFree(b);
	MatrixFree(y);
	return 0;
}
