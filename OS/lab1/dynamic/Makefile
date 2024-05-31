compilation:
	gcc -c hello.c
	gcc -fPIC -c dyn_hello.c
	gcc -shared -o libhello-dynamic.so dyn_hello.o
	gcc hello.o -L. -lhello-dynamic -o hello_dyn

moving:
	cp libhello-dynamic.so /lib/
