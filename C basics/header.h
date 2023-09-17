#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Tree {
	unsigned char symb;
	size_t freq;
	struct Tree* left;
	struct Tree* right;
}tree;

struct list {
	tree* data;
	struct list* next;
};

typedef struct stack {
	struct list* top;
}Stack;