main: main.o 
	g++ -o main main.o
generate: generate.o
	g++ -o generate generate.o
compare: compare.o
	g++ -o compare compare.o
print: print.o
	g++ -o print print.o


generate.o: generate.cpp
	g++ -c -o generate.o generate.cpp
main.o: main.cpp
	g++ -c -o main.o main.cpp
compare.o: compare.o
	g++ -c -o compare.o compare.cpp
print.o: print.cpp
	g++ -c -o print.o print.cpp

test: compare main 
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2
	./main "./test/A.dat" "./test/B.dat" "./test/C3.dat" 3
	./main "./test/A.dat" "./test/B.dat" "./test/C4.dat" 4
	./main "./test/A.dat" "./test/B.dat" "./test/C5.dat" 5
	./compare "./test/C0.dat" "./test/C.dat"
	./compare "./test/C1.dat" "./test/C.dat"
	./compare "./test/C2.dat" "./test/C.dat"
	./compare "./test/C3.dat" "./test/C.dat"
	./compare "./test/C4.dat" "./test/C.dat"
	./compare "./test/C5.dat" "./test/C.dat"
	rm -f ./test/C0.dat ./test/C1.dat ./test/C2.dat ./test/C3.dat ./test/C4.dat ./test/C5.dat

report: generate main
	rm -f dat.txt
	./generate f 50 50 A.dat 
	./generate f 50 50 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./sr
	gnuplot plotset
	mv plot.svg plot50x50.svg
	rm -f dat.txt
	./generate f 200 200 A.dat 
	./generate f 200 200 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./sr
	gnuplot plotset
	mv plot.svg plot200x200.svg
	rm -f dat.txt
	./generate f 300 300 A.dat 
	./generate f 300 300 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./sr
	gnuplot plotset
	mv plot.svg plot300x300.svg
	rm -f dat.txt
	./generate f 500 500 A.dat 
	./generate f 500 500 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./sr
	gnuplot plotset
	mv plot.svg plot500x500.svg
	rm -f dat.txt
	./generate f 1000 1000 A.dat 
	./generate f 1000 1000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	./sr
	gnuplot plotset
	mv plot.svg plot1000x1000.svg
	rm -f dat.txt A.dat B.dat C.dat



clean : 
	rm -f *.dat *.txt *.o *.svg