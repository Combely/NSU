valgrind:
	g++ main.cpp -g -o main
	valgrind --leak-check=full ./main

sanitilizer:
	g++ -fsanitize=address -fno-omit-frame-pointer -g -o mainsan main.cpp
	./mainsan 1500
