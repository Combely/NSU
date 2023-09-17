#pragma once
#include "header.h"

void encode(tree* huffTreeNode, unsigned char** huffCodes, unsigned char* currCode, size_t codeLen);

void treeSerialization(tree* node, FILE* out);

size_t getLenofOriginalFile(size_t* freqTable);

size_t getLenofEncodedFile(size_t* freqTable, unsigned char** huffCodes);

size_t getLenofTreeTopology(size_t* freqTable);

void compressing(FILE* in, FILE* out, unsigned char* inputFileName);