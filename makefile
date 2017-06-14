All: main

main: main.cpp simplex.o
	g++ main.cpp simplex.o -static -std=c++11 -Wall -Wno-unused-result -O2 -lm

simplex: simplex.cpp
	g++ -c simplex.cpp -static -std=c++11 -Wall -Wno-unused-result -O2 -lm

clean:
	rm -f *.o a.out
