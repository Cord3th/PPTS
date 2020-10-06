main: main.o
	g++ -o main main.o
generate: generate.o
	g++ -o generate generate.o
compare: compare.o
	g++ -o compare compare.o
print: print.o
	g++ -o print print.o
sr: sr.o
	g++ -o sr sr.o

generate.o: generate.cpp
	g++ -c -o generate.o generate.cpp
main.o: main.cpp
	g++ -c -o main.o main.cpp
compare.o: compare.o
	g++ -c -o compare.o compare.cpp
print.o: print.cpp
	g++ -c -o print.o print.cpp
sr.o: sr.cpp
	g++ -c -o sr.o sr.cpp

test: compare main
	for NUMBER in 0 1 2 3 4 5 ; do \
		./main "./test/A.dat" "./test/B.dat" "./test/C$$NUMBER.dat" $$NUMBER ; \
		./compare "./test/C$$NUMBER.dat" "./test/C.dat"  || exit $(.SHELLSTATUS) ; \
	done
	rm -f ./test/C0.dat ./test/C1.dat ./test/C2.dat ./test/C3.dat ./test/C4.dat ./test/C5.dat
report: generate main sr
	rm -f dat.txt
	./generate f 50 50 A.dat
	./generate f 50 50 B.dat
	for ITERATION in 0 1 2 3 4 5 6 7 8 9; do \
		for NUMBER in 0 1 2 3 4 5 ; do \
			./main A.dat B.dat C.dat $$NUMBER ; \
		done \
	done
	./sr
	gnuplot plotset
	mv plot.svg plot50x50.svg
	rm -f dat.txt
	./generate f 200 200 A.dat
	./generate f 200 200 B.dat
	for ITERATION in 0 1 2 3 4 5 6 7 8 9; do \
		for NUMBER in 0 1 2 3 4 5 ; do \
			./main A.dat B.dat C.dat $$NUMBER ; \
		done \
	done
	./sr
	gnuplot plotset
	mv plot.svg plot200x200.svg
	rm -f dat.txt
	./generate f 300 300 A.dat
	./generate f 300 300 B.dat
	for ITERATION in 0 1 2 3 4 5 6 7 8 9; do \
		for NUMBER in 0 1 2 3 4 5 ; do \
			./main A.dat B.dat C.dat $$NUMBER ; \
		done \
	done
	./sr
	gnuplot plotset
	mv plot.svg plot300x300.svg
	rm -f dat.txt
	./generate f 500 500 A.dat
	./generate f 500 500 B.dat
	for ITERATION in 0 1 2 3 4 5 6 7 8 9; do \
		for NUMBER in 0 1 2 3 4 5 ; do \
			./main A.dat B.dat C.dat $$NUMBER ; \
		done \
	done
	./sr
	gnuplot plotset
	mv plot.svg plot500x500.svg
	rm -f dat.txt
	./generate f 1000 1000 A.dat
	./generate f 1000 1000 B.dat
	for ITERATION in 0 1 2 3 4 5 6 7 8 9; do \
		for NUMBER in 0 1 2 3 4 5 ; do \
			./main A.dat B.dat C.dat $$NUMBER ; \
		done \
	done
	./sr
	gnuplot plotset
	mv plot.svg plot1000x1000.svg
	rm -f dat.txt A.dat B.dat C.dat

clean:
	rm -f *.dat *.txt *.o *.svg
