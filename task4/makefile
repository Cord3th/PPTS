main: main.cpp
	mpicxx -O3 -o main main.cpp
generate: generate.cpp
	g++ -O3 -o generate generate.cpp
print: print.cpp
	g++ -O3 -o print print.cpp
random: random.cpp
	g++ -O3 -o random random.cpp


genmatr: generate
	./generate f 10 10 A.dat
	./generate f 10 1 b.dat

test:generate print main
	./generate f 10 10 A.dat
	./generate f 10 1 b.dat
	./print A.dat A.txt
	./print b.dat b.txt
	mpirun -n 5 main A.dat b.dat c.dat

clean : 
	rm -f *.dat *.txt  *.svg print generate main random *.map