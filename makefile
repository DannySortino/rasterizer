make: main.cpp
	g++ -o main main.cpp bmpLibrary/EasyBMP.cpp -std=c++11 -Wall

run: main.cpp
	rm -f main; g++ -o main main.cpp bmpLibrary/EasyBMP.cpp; ./main
