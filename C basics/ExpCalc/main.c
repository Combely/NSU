#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>

double exp(double x, double accuracy) {
	double sum = 0, locsum = 1;
	int n = 1;
	while (locsum > accuracy) {
		sum += locsum;
		locsum *= x / n;
		n++;
	}
	return sum;
}

int main() {
	FILE* in = fopen("input.txt", "r"), * out = fopen("output.txt", "w");
	double E = 0.000000000001, currx, sum = 0;
	int n = 1, N;
	fscanf(in, "%d", &N);
	for (int i = 0; i < N; i++) {
		fscanf(in, "%lf", &currx);
		//в условии сказано не подключать math.h, fabs из math.h, а простой abs не подойдет, т.к. он для int и округлит вещественные
		sum = exp(currx > 0.0 ? currx : -currx, E);
		sum = currx > 0.0 ? sum : 1 / sum;
		fprintf(out, "%0.15g\n", sum);
	}
	return 0;
}