#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>
#include <sys/errno.h>

#define PATH_LEN 512

const char basic_path [PATH_LEN] = "/home/combely/NSUGit/OS/lab3/";
const char SRCdir_name [PATH_LEN] = "source";

char* reverseStr(char* str);
void reverseFile(char *name, char *workPath, char *newPath);
void reverseDir(char *name, char *workPath, char *newPath);
void createRevDir(char* name, char* path);

int main() {
    char* SRCpath = (char*)malloc(PATH_LEN * sizeof(char));
    char* DSTpath = (char*)malloc(PATH_LEN * sizeof(char));
    strcpy(SRCpath, basic_path);
    strcpy(DSTpath, basic_path);
	reverseDir(SRCdir_name, SRCpath, DSTpath);
	return 0;
}

void reverseDir(char *dir_name, char *SRCpath, char *DSTpath) {
	printf("Mirrorcopying %s from %s\ninto %s\n", dir_name, SRCpath, DSTpath);
	DIR *current = opendir(strcat(SRCpath, dir_name)); //will get "/home/combely/NSUGit/OS/lab3/source"
	struct dirent *entry;
	char *curPath = (char*)malloc(PATH_LEN * sizeof(char));
	strcpy(curPath, SRCpath);
	strcat(curPath,"/");
	createDir(dir_name, DSTpath);

	while ((entry = readdir(current)) != NULL) {
		if (entry->d_type == DT_REG)
			reverseFile(entry->d_name, curPath, DSTpath);
		else if (entry->d_type == DT_DIR)
			if (strcmp(entry->d_name, "..") && strcmp(entry->d_name, "."))
				reverseDir(entry->d_name, curPath, DSTpath);
	}
	printf("Mirrorcopying %s finished\n\n", dir_name);
	closedir(current);
	free(entry);
}

char* reverseSTR(char* basic) {
	int len = strlen(basic);
	char* reversed = (char*)malloc(len * sizeof(char));
	for (int i = 0; i < len; i++)
		reversed[i] = basic[len - i - 1];
	return reversed;
}

void reverseFile(char *name, char *SRCpath, char *DSTpath){
	printf("Rewriting file %s from %s\ninto %s\n", name, SRCpath, DSTpath);
	char c;
	char *buf = reverseSTR(name);
	char *tmp = (char*)malloc(PATH_LEN * sizeof(char));
	strcpy(tmp, SRCpath);
	
	FILE *oldFile;
	if (oldFile = fopen(strcat(tmp, name), "rb")) {
		printf("Error during opening of file %s ", tmp);
		if (errno == EACCES) printf("(access denied)\n");
		if (errno == ENOENT) printf("(incorrect path)\n");
		exit(1);
	}

	strcpy(tmp, DSTpath);
	FILE *newFile = fopen(strcat(tmp, buf), "w");
	int i = -1;

	while (fseek(oldFile, i, SEEK_END) != -1) {
		c = fgetc(oldFile);
		// printf("%c", c);
		fputc(c, newFile);
		--i;
	}
	if (fclose(oldFile)) {
		printf("Error during closing of file %s", tmp);
		if (errno == EBADF) printf("(bad file descriptor)\n");
		exit(1);
	}
	fclose(newFile);
	free(buf);
	free(tmp);
	printf("Rewriting file %s succeeded\n\n", name, SRCpath, DSTpath);
}

void createDir(char* name, char* path)
{
	char* buf = reverseStr(name);
	strcat(path, buf);
	mkdir(path, 0777);
	strcat(path, "/");
	free(buf);
}
