#include <stdio.h>

int main() {
	FILE *pfile = fopen("hello_out.txt", "w");
	fprintf(pfile, "Hello world!\n");
	fclose(pfile);
	return 0;
}
