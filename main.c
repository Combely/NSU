#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc.h"

void inf2post(char* outstr, FILE* in);

void stackcalc(char* inputstr);

void UI(FILE **in) {
	printf("Description of input values:");
	printf("\nAn expression with a unary minus must be enclosed in parentheses on both sides if there is an operation sign in front of it");
	printf("\nUse command -s to set input file default \"input.txt\": ");
	char str[3] = { 0 }, *pstr;
	pstr = gets(str);
	if (pstr[0] == '-' && pstr[1] == 's')
		in = fopen("input.txt", "r");
	else {
		printf("poshel naxyi");
	}
}

int main(int argc, char **argv) {
	printf("ARGC = %d\n", argc);
	FILE* in;
	char str[500001] = { 0 };
	if (argc > 1) {
		in = fopen(argv[1], "r");
	}
	else {
		UI(&in);
	}
	inf2post(str, in);
	stackcalc(str);
	return 0;
}