#include <stdio.h>
#include "static_hello.h"

void print_hello() {
    printf("hello world!\n");
}

int main () {
    print_hello();
    hello_from_static_lib();
    return 0;
}
