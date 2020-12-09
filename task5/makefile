main: main.cpp
	mpicxx -std=c++11 -O3 -o main main.cpp
generate: generate.cpp
	g++ -std=c++11 -O3 -o generate generate.cpp
print: print.cpp
	g++ -std=c++11 -O3 -o print print.cpp
compare: compare.cpp
	g++ -std=c++11 -O3 -o compare compare.cpp


gen: generate
	./generate d 10 10 A.dat
	./generate d 10 10 B.dat
printmat: print
	./print A.dat A.txt
	./print B.dat B.txt

test: main compare
	mpirun -np 8 main "./test/A.dat" "./test/B.dat" "./test/C0.dat"
	mpirun -np 27 main "./test/A.dat" "./test/B.dat" "./test/C1.dat"
	./compare "./test/C0.dat" "./test/C.dat"
	./compare "./test/C1.dat" "./test/C.dat"
	rm "./test/C0.dat"
	rm "./test/C1.dat"
clean:
	rm -f *.dat *.txt  *.svg print generate main random compare *.map
