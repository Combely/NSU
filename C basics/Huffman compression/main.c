#include "header.h"
#include "auxillaryfuncs.h"
#include "compressing.h"
#include "decompressing.h"

void UI(FILE** file, FILE** arc, int argc, char** argv) {
	size_t iterator = 1;
	//Âûâåñòè ñïðàâî÷íóþ èíôîðìàöèþ
	char hcnd[3] = "-h";
	//Âûâåñòè èíôîðìàöèþ èç çàãîëîâêà àðõèâà ñ èìåíåì arc
	char icnd[3] = "-i";
	//Èçâëå÷ü ôàéë ñ íàçâàíèåì file èç àðõèâà arc
	char xcnd[3] = "-x";
	//Ïîìåñòèòü ôàéë file â àðõèâ arc
	char acnd[3] = "-a";
	if (argc == 1) {
		printf("Enter the command -h for reference information");
	}
	else {
		while (1) {
			if (!strcmp(argv[iterator], hcnd)) {
				printf("Commands:\n");
				printf("\t -i [name of archive] \t The program will show information about compressed file such as Huffman tree and compressing ratio.\n");
				printf("\t -x [name of file] [name of arcive] \t First argument - name of file to extract from arcive with name of Second argument. \n");
				printf("\t -a [name of file] [name of arcive] \t First argument - name of file to compress. Second argument - name of arcive to write to compressed file to.\n");
				printf("Tip to correct display of decompressed files:\n");
				printf("Enter the file extension after the file name (commands -a, -x)\n");
				iterator++;
			}
			else if (!strcmp(argv[iterator], icnd)) {
				*arc = fopen(argv[iterator + 1], "rb");
				if (*arc == NULL)
					printf("Archive \"%s\" open error : no such file in the directory", argv[iterator + 1]);
				else
					arcInfo(*arc);
				iterator += 2;
			}
			else if (!strcmp(argv[iterator], xcnd)) {
				*file = fopen(argv[iterator + 1], "wb"); 
				*arc = fopen(argv[iterator + 2], "rb");
				if (*arc == NULL)
					printf("Archive \"%s\" open error : no such file in the directory", argv[iterator + 2]);
				else {
					char check;
					if (fread(&check, 1, 1, *arc) == 1) {
						rewind(*arc);
						decompressing(*arc, *file, argv[iterator + 1]);
						printf("Decompression passed succesfully\n");
					}
					else {
						printf("Decompression can not be percised: archive \"%s\" is empty", argv[iterator + 2]);
					}
				}
				iterator += 3;
			}
			else if (!strcmp(argv[iterator], acnd)) {
				*file = fopen(argv[iterator + 1], "rb");
				*arc = fopen(argv[iterator + 2], "wb");
				if (*file == NULL)
					printf("File \"%s\" open error: no such file in the directory", argv[iterator + 1]);
				else {
					char check;
					if (fread(&check, 1, 1, *file) == 1) {
						rewind(*file);
						compressing(*file, *arc, argv[iterator + 1]);
						printf("Compression passed succesfully\n");
					}
					else {
						printf("Compression can not be percised: file \"%s\" is empty", argv[iterator + 1]);
					}
				}
				iterator += 3;
			}
			else { //Âî èçáåæàíèå çàöèêëèâàíèÿ â ñëó÷àå íåâåðíî ââåäåííûõ êîìàíä
				iterator++;
				printf("Wrong input. Enter the command -h for reference information.");
			}

			if (iterator >= argc)
				break;
		}
	}
}

int main(int argc, char** argv) {
	FILE* in = NULL, * arc = NULL;
	UI(&in, &arc, argc, argv);
	return 0;
}

