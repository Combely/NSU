#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc.h"

Stack* create() {
	Stack* S = (Stack*)malloc(sizeof(Stack));
	S->top = NULL;
	return S;
}

void push(Stack* S, double a) {
	struct list* p = (struct list*)malloc(sizeof(struct list));
	p->data = a;
	p->next = S->top;
	S->top = p;
}

int empty(Stack* S) {
	return (S->top == NULL);
}

double top(Stack* S) {
	return (S->top->data);
}

double pop(Stack* S) {
	struct list* t = S->top;
	double a = t->data;
	S->top = t->next;
	free(t);
	return a;
}

int priority(char a, char b) {
	int priority[2] = { 0 }, i;
	char curr;
	while (priority[0] == 0 || priority[1] == 0) {
		if (priority[0] == 0) {
			curr = a;
			i = 0;
		}
		else {
			curr = b;
			i = 1;
		}
		if (curr == '(')
			priority[i] = 1;
		else if (curr == ')')
			priority[i] = 2;
		else if (curr == '=')
			priority[i] = 3;
		else if (curr == '+' || curr == '-')
			priority[i] = 4;
		else if (curr == '*' || curr == '/')
			priority[i] = 5;
	}
	return (priority[0] >= priority[1]);
}

void strapend(char* str, char symb, int len) {
	str[len] = symb;
	str[len + 1] = '\0';
}

void inf2post(char* outstr, FILE *in) {
	Stack* S = create();
	int pty, mainind = 0, ubminchecker = 1, outlen = 0; //åñëè ôëàã 1, òî ìèíóñ óíàðíûé, åñëè ôëàã 0 - áèíàðíûé
	char* num, * pinstr, * instr;
	double trash;
	instr = (char*)malloc(500020 * sizeof(char));
	pinstr = fgets(instr, 500020, in);
	while (pinstr[mainind] != '\0') {
		if (pinstr[mainind] == ' ') {
			mainind++;
			continue;
		}
		else if (pinstr[mainind] >= '0' && pinstr[mainind] <= '9') {
			strapend(outstr, pinstr[mainind], outlen);
			outlen++;
			if (pinstr[mainind + 1] >= '0' && pinstr[mainind + 1] <= '9') {
				mainind++;
				continue;
			}
			else {
				strapend(outstr, ' ', outlen);
				outlen++;
			}
			ubminchecker = 0;
		}
		else if (pinstr[mainind] == '(') {
			push(S, pinstr[mainind]);
			ubminchecker = 1;
		}
		else if (pinstr[mainind] == ')') {
			while ((int)top(S) != '(') {
				strapend(outstr, (int)pop(S), outlen);
				outlen++;
				strapend(outstr, ' ', outlen);
				outlen++;
			}
			trash = pop(S);
		}
		else {
			if (ubminchecker) {
				strapend(outstr, '0', outlen);
				outlen++;
				strapend(outstr, ' ', outlen);
				outlen++;
				ubminchecker = 0;
				continue;
			}
			else {
				if (!empty(S))
					pty = priority((int)top(S), pinstr[mainind]);
				while (!empty(S) && pty) {
					strapend(outstr, (int)pop(S), outlen);
					outlen++;
					strapend(outstr, ' ', outlen);
					outlen++;
					if (!empty(S))
						pty = priority((int)top(S), pinstr[mainind]);
				}
				push(S, pinstr[mainind]);
			}
		}
		mainind++;
	}
	while (!empty(S)) {
		strapend(outstr, (int)pop(S), outlen);
		outlen++;
		strapend(outstr, ' ', outlen);
	}
	free(S);
	free(instr);
}

void stackcalc(char* inputstr) {
	Stack* S = create();
	char* num, * trash;
	int mainind = 0, len = 0;
	double rightop, leftop, ans = 0, token;
	num = (char*)malloc(500020 * sizeof(char));
	num[0] = 0;
	while (inputstr[mainind] != '\0') {
		if (inputstr[mainind] == ' ') {
			mainind++;
			continue;
		}
		else if (inputstr[mainind] == '+' || inputstr[mainind] == '-' || inputstr[mainind] == '*' || inputstr[mainind] == '/') {
			rightop = pop(S);
			leftop = pop(S);
			if (inputstr[mainind] == '+')
				ans = leftop + rightop;
			else if (inputstr[mainind] == '-')
				ans = leftop - rightop;
			else if (inputstr[mainind] == '*')
				ans = leftop * rightop;
			else if (inputstr[mainind] == '/') {
				ans = leftop / rightop;
			}
			push(S, ans);
		}
		else {
			strapend(num, inputstr[mainind], len);
			len++;
			if (inputstr[mainind + 1] >= '0' && inputstr[mainind + 1] <= '9') {
				mainind++;
				continue;
			}
			else {
				token = strtod(num, &trash);
				push(S, token);
				num[0] = 0; //î÷èùàåì num
				len = 0;
			}
		}
		mainind++;
	}
	printf("%0.20lf", pop(S));
	free(S);
	free(num);
}
