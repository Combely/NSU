#pragma once
#include "header.h"

tree* treeDeserialization(FILE* in, size_t treeTopologyLen);

void decompressing(FILE* arc, FILE* decodedFile, unsigned char* nameofEncodedFile);
