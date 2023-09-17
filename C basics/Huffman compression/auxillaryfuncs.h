#pragma once
#include "header.h"

void arcInfo(FILE* arc);

void freeTree(tree* node);

tree* createNode(unsigned char symb, size_t freq, tree* leftChild, tree* rightChild);

Stack* createStack();

void push(Stack* S, tree* a);

size_t empty(Stack* S);

tree* top(Stack* S);

tree* pop(Stack* S);

void siftup(tree** mas, size_t ind);

void insert(tree** heap, size_t* size, tree* curr);

size_t parentheap(size_t ind);

size_t leftheap(size_t ind);

size_t rightheap(size_t ind);

void siftdown(tree** mas, size_t ind, size_t size);

tree* extractmin(tree** heap, size_t* size);

void insert(tree** heap, size_t* size, tree* curr);

void siftup(tree** mas, size_t ind);