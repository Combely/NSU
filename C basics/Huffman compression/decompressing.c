#include "header.h"
#include "auxillaryfuncs.h"

tree* treeDeserialization(FILE* in, size_t treeTopologyLen) {
	Stack* stack = createStack();
	tree* root;
	size_t iterator = 0;
	unsigned char cur;
	while (iterator < treeTopologyLen) {
		fread(&cur, sizeof(unsigned char), 1, in);
		if (cur == '1') {
			tree* newNode;
			fread(&cur, sizeof(unsigned char), 1, in);
			//Из-за особенности построения дерева в compressing на данном этапе величина node->freq лишняя, поэтому
			//инициализируем ее фиктивными нулями.
			newNode = createNode(cur, 0, NULL, NULL);
			push(stack, newNode);
			iterator += 2;
		}
		else if (cur == '0') {
			tree* rightsym = pop(stack);
			tree* leftsym = pop(stack);
			tree* newNode = createNode(0, 0, leftsym, rightsym);
			push(stack, newNode);
			iterator++;
		}
	}
	root = pop(stack);
	free(stack);
	return root;
}

void decompressing(FILE* arc, FILE* decodedFile, unsigned char* nameofEncodedFile) {
	tree* root, * curnode;
	size_t treeLen, nameofFileLen, encodedfileLen, curbit;
	unsigned char buf;
	fseek(arc, sizeof(size_t), SEEK_CUR);
	fread(&encodedfileLen, sizeof(size_t), 1, arc);
	fread(&nameofFileLen, sizeof(size_t), 1, arc);
	fseek(arc, sizeof(unsigned char) * nameofFileLen, SEEK_CUR);
	fread(&treeLen, sizeof(size_t), 1, arc);
	root = treeDeserialization(arc, treeLen);
	curnode = root;
	//Запись декодированного буфера 
	while (encodedfileLen) {
		fread(&buf, sizeof(unsigned char), 1, arc);
		for (int bufbitsInuse = 7; bufbitsInuse >= 0; bufbitsInuse--) {
			if (!encodedfileLen)
				break;
			curbit = 1 & (buf >> bufbitsInuse);
			//Для случая, когда файл состоит из подряд написанного единственного символа
			if (!curnode->left && !curnode->right) {
				fwrite(&curnode->symb, sizeof(unsigned char), 1, decodedFile);
				curnode = root;
			}
			else {
				if (curbit == 1)
					curnode = curnode->right;
				else if (curbit == 0)
					curnode = curnode->left;
				//Если после перемещения ниже по дереву попали в лист, то сразу выводим соответсвующий символ
				if (!curnode->left && !curnode->right) {
					fwrite(&curnode->symb, sizeof(unsigned char), 1, decodedFile);
					curnode = root;
				}
			}
			encodedfileLen--;
		}
	}
}