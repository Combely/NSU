#include <stdio.h>
#include "dyn_hello.h"

void print_hello() {
    printf("hello world!\n");
}

int main () {
    print_hello();
    hello_from_dynamic_lib();
    return 0;
}
