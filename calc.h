
#ifndef CALC_H
#define CALC_H

struct list {
	double data;
	struct list* next;
};

typedef struct stack {
	struct list* top;
} Stack;

Stack* create();

void push(Stack* S, double a);

int empty(Stack* S);

double top(Stack* S);

double pop(Stack* S);

int priority(char a, char b);

void strapend(char* str, char symb, int len);

void inf2post(char* outstr, FILE *in);

void stackcalc(char* inputstr);

#endif
