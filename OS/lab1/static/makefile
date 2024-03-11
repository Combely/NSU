static:
	gcc -c hello.c
	gcc -c static_hello.c
	ar rc libhello-static.a static_hello.o
	ranlib libhello-static.a
	gcc hello.o -L. -lhello-static -o hello_static
