#include "header.h"
#include "auxillaryfuncs.h"


void encode(tree* huffTreeNode, unsigned char** huffCodes, unsigned char* currCode, size_t codeLen) {
	if (!huffTreeNode)
		return;
	if (!huffTreeNode->left && !huffTreeNode->right) {
		if (!codeLen) { //��� ������, ����� ���� �������� 1 ���������� ������
			huffCodes[huffTreeNode->symb] = "0";
		}
		else 
			huffCodes[huffTreeNode->symb] = currCode;
		return;
	}
	unsigned char* codeForLeft = (unsigned char*)malloc((codeLen + 2) * sizeof(unsigned char));
	unsigned char* codeForRight = (unsigned char*)malloc((codeLen + 2) * sizeof(unsigned char));
	strcpy(codeForLeft, currCode);
	strcpy(codeForRight, currCode);
	free(currCode);
	codeForLeft[codeLen] = '0';
	codeForRight[codeLen] = '1';
	codeForLeft[codeLen + 1] = '\0';
	codeForRight[codeLen + 1] = '\0';
	encode(huffTreeNode->left, huffCodes, codeForLeft, codeLen + 1);
	encode(huffTreeNode->right, huffCodes, codeForRight, codeLen + 1);
}

void treeSerialization(tree* node, FILE* out) {
	unsigned char zerobuf = '0';
	if (!node->left && !node->right) {
		unsigned char buf = '1';
		fwrite(&buf, sizeof(unsigned char), 1, out);
		fwrite(&(node->symb), sizeof(unsigned char), 1, out);
		return;
	}
	treeSerialization(node->left, out);
	treeSerialization(node->right, out);
	fwrite(&zerobuf, sizeof(unsigned char), 1, out);
}

size_t getLenofOriginalFile(size_t* freqTable) {
	size_t len = 0;
	for (size_t i = 0; i < 256; i++) {
		if (freqTable[i])
			len += freqTable[i];
	}
	return len;
}

size_t getLenofEncodedFile(size_t* freqTable, unsigned char** huffCodes) {
	size_t len = 0;
	for (size_t i = 0; i < 256; i++) {
		if (freqTable[i])
			len += (size_t)strlen(huffCodes[i]) * freqTable[i];
	}
	return len;
}

//���������� ���������� ����, ������� ������ ������
size_t getLenofTreeTopology(size_t* freqTable) {
	size_t len = 0, temp;
	for (size_t i = 0; i < 256; i++) {
		if (freqTable[i]) {
			len++;
		}
	}
	temp = len;
	temp--;
	len *= 2;
	len += temp;
	return len;
}

void compressing(FILE* in, FILE* out, unsigned char* inputFileName) {
	unsigned char symb, ** huffCodes, * currCode, buf = 0;
	size_t* freqtable, heapsize = 0, currCodeLen = 0, bitcter = 0;
	size_t originalFLen, encodedFLen, treeLen, nameLen;
	tree** heap, * root;
	huffCodes = (unsigned char**)malloc(256 * sizeof(unsigned char*));
	currCode = (unsigned char*)malloc(sizeof(unsigned char));
	currCode[0] = '\0';
	freqtable = (unsigned int*)malloc(256 * sizeof(unsigned int));
	heap = (tree**)malloc(256 * sizeof(tree*));
	for (size_t i = 0; i < 256; i++) {
		huffCodes[i] = NULL;
		freqtable[i] = 0;
	}
	//������� ������ ��������� �������� � �����
	while (fread(&symb, sizeof(unsigned char), 1, in) == 1) {
		freqtable[(int)symb]++;
	}
	//���������� ������ �������� 
	for (size_t i = 0; i < 256; i++) {
		//������� ��� ������� � ���������� ��������� � ����, ����� �� ��������� ����� ������ ������ � �� ����������� ���
		if (freqtable[i]) {
			tree* newNode = createNode((unsigned char)i, freqtable[i], NULL, NULL);
			insert(heap, &heapsize, newNode);
		}
	}
	while (heapsize != 1) {
		//���� � ���� �� ��������� ������ ������� ������, ��������� ��� �������� � ������������ ���������,
		//���������� �� ��� ����� ���� � ��������, ������ ����� ������ ��� �������
		tree* left = extractmin(heap, &heapsize);
		tree* right = extractmin(heap, &heapsize);
		tree* newNode = createNode('\0', (left->freq + right->freq), left, right);
		insert(heap, &heapsize, newNode);
	}
	//��������� �� ���� ��������� �� ������� ������ ��������
	root = extractmin(heap, &heapsize);
	//���������� �� ������ ����� �������� � ������� ��������� �������
	encode(root, huffCodes, currCode, currCodeLen);
	//������� ����� ��������� � ������� �����, ����� ������ ������
	originalFLen = getLenofOriginalFile(freqtable);
	encodedFLen = getLenofEncodedFile(freqtable, huffCodes);
	treeLen = getLenofTreeTopology(freqtable);
	//������ ����� ����� � ��������� ����� ������ ��� ��������� �� ���������� �������������� ������
	nameLen = strlen(inputFileName);
	//������ ���������� � ������
	fwrite(&originalFLen, sizeof(size_t), 1, out);
	fwrite(&encodedFLen, sizeof(size_t), 1, out);
	fwrite(&nameLen, sizeof(size_t), 1, out);
	fwrite(inputFileName, sizeof(unsigned char), nameLen, out);
	fwrite(&treeLen, sizeof(size_t), 1, out);
	//������ ���������� � ������ ��� ������������� ������ 
	treeSerialization(root, out);

	/*����� �������, � ��������� ����������:
	* size_t ����� ��������� ����� (� ������)
	* size_t ����� ��������������� ����� (� �����, � ����� ��������� ������ �� ����� � ��������� �����)
	* size_t ����� �������� ����� (� ������)
	* unsigned char * (����� ��������) �������� �����
	* size_t ����� ������ ��������� ������ ��������
	* unsigned char * (����� ������ ���������) ������ ���������� ������
	������ �������������� � ������� �������� �������� ������: ����������� ����������� ����� ���������.
	���� �������� ����, ������������ 1 � ��������������� ����� ������.
	���� �������� ����, �� ���������� ������, ������������ 0.
	��� ������ ������������ ����� (� �� ����, ��� ��� �������� � ������ ����������� �����). ����� ��� �������������
	������ ����� �������� ����: 2 * (���������� �������) + (���������� ����� � ���������).
	*/

	rewind(in);
	//������ ��������������� ������
	while (fread(&symb, sizeof(unsigned char), 1, in) == 1) {
		size_t len = strlen(huffCodes[symb]);
		for (size_t i = 0; i < len; i++) {
			//1 � 0 � ������ ���� �������� �������� � ���� unsigned char, ��� �������� � int �������� ������� ���
			//����� �������, ������� ���������� �������� �� ���� '0', �.�. 49 �� ASCII
			buf |= (huffCodes[symb][i] - '0') << (7 - bitcter);
			bitcter++;
			if (bitcter == 8) {
				bitcter = 0;
				fwrite(&buf, sizeof(unsigned char), 1, out);
				buf = 0;
			}
		}
	}
	//������ ������ � �������
	if (bitcter) 
		fwrite(&buf, sizeof(unsigned char), 1, out);
	free(heap);
	free(freqtable);
	free(huffCodes);
	freeTree(root);
}
