#include "header.h"


void arcInfo(FILE* arc) {
	size_t originalFileLen, encodedfileLen, treeLen, nameofFileLen, EFLinBytes;
	double comprRatio; 
	unsigned char* treeTopology, *nameofFile;
	fread(&originalFileLen, sizeof(size_t), 1, arc);
	fread(&encodedfileLen, sizeof(size_t), 1, arc);
	fread(&nameofFileLen, sizeof(size_t), 1, arc);
	nameofFile = (unsigned char*)malloc((nameofFileLen + 1) * sizeof(unsigned char));
	fread(nameofFile, sizeof(unsigned char), nameofFileLen, arc); 
	nameofFile[nameofFileLen] = '\0';
	fread(&treeLen, sizeof(size_t), 1, arc);
	treeTopology = (unsigned char*)malloc((treeLen + 1) * sizeof(unsigned char));
	fread(treeTopology, sizeof(unsigned char), treeLen, arc);
	treeTopology[treeLen] = '\0';
	comprRatio = ((double)(originalFileLen * 8) / (double)encodedfileLen) * 100;
	EFLinBytes = (size_t)(encodedfileLen / 8);
	if ((encodedfileLen % 8))
		EFLinBytes++;
	printf("Size of original file \"%s\" is %u byte(s)\n", nameofFile, originalFileLen);
	printf("Size of compressed file \"%s\" is %u byte(s) (%u bits)\n", nameofFile, EFLinBytes, encodedfileLen);
	printf("Size of Huffman tree %u byte(s) \n", treeLen);
	//printf("Huffman tree for file \"%s\" : %s \n", nameofFile, treeTopology);
	printf("Compressing ratio of file \"%s\" is %.5lf%% ", nameofFile, comprRatio);
	printf("where compressing ratio = (original file length / compressed file length) * 100%%");
}

void freeTree(tree* node) {
	if (!node->left && !node->right) {
		free(node);
		return;
	}
	//В дереве Хаффмана каждый узел имеет или 0 или 2 сыновей, поэтому проверка на существование левого и правого
	//сыновей не нужна, поскольку узел на этом этапе уже прошел проверку на одно из условий.
	freeTree(node->left);
	freeTree(node->right);
	free(node);
}

tree* createNode(unsigned char symb, size_t freq, tree* leftChild, tree* rightChild) {
	tree* new = (tree*)malloc(sizeof(tree));
	new->symb = symb;
	new->freq = freq;
	new->left = leftChild;
	new->right = rightChild;
	return new;
}

Stack* createStack() {
	Stack* S = (Stack*)malloc(sizeof(Stack));
	S->top = NULL;
	return S;
}

void push(Stack* S, tree* a) {
	struct list* p = (struct list*)malloc(sizeof(struct list));
	p->data = a;
	p->next = S->top;
	S->top = p;
}

size_t empty(Stack* S) {
	return (S->top == NULL);
}

tree* top(Stack* S) {
	return (S->top->data);
}

tree* pop(Stack* S) {
	struct list* t = S->top;
	tree* a = t->data;
	S->top = t->next;
	free(t);
	return a;
}

void siftup(tree** mas, size_t ind);
void insert(tree** heap, size_t* size, tree* curr);

size_t parentheap(size_t ind) {
	return (ind - 1) / 2;
}

size_t leftheap(size_t ind) {
	return ind * 2 + 1;
}

size_t rightheap(size_t ind) {
	return ind * 2 + 2;
}

void siftdown(tree** mas, size_t ind, size_t size) {
	size_t l = leftheap(ind), r = rightheap(ind), min = ind;
	tree* t;
	if (l < size && mas[l]->freq < mas[min]->freq)
		min = l;
	if (r < size && mas[r]->freq < mas[min]->freq)
		min = r;
	if (mas[min]->freq != mas[ind]->freq) {
		t = mas[ind];
		mas[ind] = mas[min];
		mas[min] = t;
		siftdown(mas, min, size);
	}
}

tree* extractmin(tree** heap, size_t* size) {
	tree* t;
	(*size)--;
	t = heap[0];
	heap[0] = heap[*size];
	heap[*size] = t;
	siftdown(heap, 0, *size);
	return heap[*size];
}

void insert(tree** heap, size_t* size, tree* curr) {
	(*size)++;
	heap[*size - 1] = curr;
	siftup(heap, *size - 1);
}

void siftup(tree** mas, size_t ind) {
	size_t par;
	if (!ind)
		return;
	par = parentheap(ind);
	if (mas[ind]->freq < mas[par]->freq) {
		tree* t = mas[ind];
		mas[ind] = mas[par];
		mas[par] = t;
		siftup(mas, par);
	}
}
