#include <stdio.h>
#include <sys/syscall.h> 
#include <unistd.h>

ssize_t my_write_wrapper(int file_descriptor, char* buf, size_t buf_size) {
	return syscall(SYS_write, file_descriptor, buf, buf_size);
}

int main() {
	my_write_wrapper(1, "Hello world!\n", 13);
	return 0;
}
