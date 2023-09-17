#include "header.h"
#include "auxillaryfuncs.h"


void encode(tree* huffTreeNode, unsigned char** huffCodes, unsigned char* currCode, size_t codeLen) {
	if (!huffTreeNode)
		return;
	if (!huffTreeNode->left && !huffTreeNode->right) {
		if (!codeLen) { //для случая, когда файл содержит 1 уникальный символ
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

//Возвращает количество байт, которые займет дерево
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
	//Подсчет частот вхождения символов в файле
	while (fread(&symb, sizeof(unsigned char), 1, in) == 1) {
		freqtable[(int)symb]++;
	}
	//Построение дерева Хаффмана 
	for (size_t i = 0; i < 256; i++) {
		//Заносим все символы с ненулевыми частотами в кучу, чтобы не создавать новый массив ссылок и не сортировать его
		if (freqtable[i]) {
			tree* newNode = createNode((unsigned char)i, freqtable[i], NULL, NULL);
			insert(heap, &heapsize, newNode);
		}
	}
	while (heapsize != 1) {
		//Пока в куче не останется только вершина дерева, извлекаем два элемента с минимальными частотами,
		//составляем из них новый узел с частотой, равной сумме частот его сыновей
		tree* left = extractmin(heap, &heapsize);
		tree* right = extractmin(heap, &heapsize);
		tree* newNode = createNode('\0', (left->freq + right->freq), left, right);
		insert(heap, &heapsize, newNode);
	}
	//Извлекаем из кучи указатель на вершину дерева Хаффмана
	root = extractmin(heap, &heapsize);
	//Вычисление по дереву кодов Хаффмана с помощью отдельной функции
	encode(root, huffCodes, currCode, currCodeLen);
	//Считаем длину исходного и сжатого файла, длину записи дерева
	originalFLen = getLenofOriginalFile(freqtable);
	encodedFLen = getLenofEncodedFile(freqtable, huffCodes);
	treeLen = getLenofTreeTopology(freqtable);
	//Запись имени файла в заголовок части архива для навигации по нескольким закодированным файлам
	nameLen = strlen(inputFileName);
	//Запись информации о длинах
	fwrite(&originalFLen, sizeof(size_t), 1, out);
	fwrite(&encodedFLen, sizeof(size_t), 1, out);
	fwrite(&nameLen, sizeof(size_t), 1, out);
	fwrite(inputFileName, sizeof(unsigned char), nameLen, out);
	fwrite(&treeLen, sizeof(size_t), 1, out);
	//Запись информации о дереве для декодирования буфера 
	treeSerialization(root, out);

	/*Таким образом, в заголовке содержится:
	* size_t Длина исходного файла (в байтах)
	* size_t Длина закодированного файла (в битах, с целью выявления хвоста из битов в последнем байте)
	* size_t Длина названия файла (в байтах)
	* unsigned char * (длина названия) Название файла
	* size_t Длина записи топологии дерева Хаффмана
	* unsigned char * (длина записи топологии) Запись топологиии дерева
	Дерево представляется с помощью обратной польской записи: совершается постфиксный обход структуры.
	Если встречен лист, записывается 1 и соответствующий листу символ.
	Если встречен узел, не являющийся листом, записывается 0.
	Для записи используются байты (а не биты, как это делается в случае кодирования файла). Всего для представления
	дерева будет записано байт: 2 * (количество листьев) + (количество узлов с сыновьями).
	*/

	rewind(in);
	//Запись закодированного буфера
	while (fread(&symb, sizeof(unsigned char), 1, in) == 1) {
		size_t len = strlen(huffCodes[symb]);
		for (size_t i = 0; i < len; i++) {
			//1 и 0 в строке кода Хаффмана записаны в виде unsigned char, при переводе в int оператор возьмет код
			//этого символа, поэтому изначально вычитаем из кода '0', т.е. 49 из ASCII
			buf |= (huffCodes[symb][i] - '0') << (7 - bitcter);
			bitcter++;
			if (bitcter == 8) {
				bitcter = 0;
				fwrite(&buf, sizeof(unsigned char), 1, out);
				buf = 0;
			}
		}
	}
	//Запись буфера с хвостом
	if (bitcter) 
		fwrite(&buf, sizeof(unsigned char), 1, out);
	free(heap);
	free(freqtable);
	free(huffCodes);
	freeTree(root);
}
