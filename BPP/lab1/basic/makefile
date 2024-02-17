valgrind:
	g++ grad.cpp -g -o grad
	valgrind --leak-check=full ./grad

sanitilizer:
	g++ -fsanitize=address -fno-omit-frame-pointer -g -o gradsan grad.cpp
	./gradsan 1000
