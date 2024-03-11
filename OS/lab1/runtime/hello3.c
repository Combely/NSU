#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "dyn_runtime_hello.h"

void print_hello() {
	printf("hello world!\n");
}

int main () {
    	void *handle;
    	char *error;
	void (*hellofunc_addr)(void);
    	handle = dlopen ("/home/combely/OSlab1/libhello-dynruntime.so", RTLD_LAZY);
	if (!handle) {
		fprintf(stderr, "%s\n", dlerror());
		exit(EXIT_FAILURE);
	}
	dlerror();//clear existing error
	*(void **) (&hellofunc_addr) = dlsym(handle, "hello_from_dyn_runtime_lib");
	error = dlerror();
	if (error != NULL) {
		fprintf(stderr, "%s\n", error);
		exit(EXIT_FAILURE);
	}
    	print_hello();
    	(*hellofunc_addr)();
    	dlclose(handle);
	exit(EXIT_SUCCESS);
}
