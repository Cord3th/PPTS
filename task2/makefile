.PHONY: test plot

main: main.cpp
	g++ -std=c++11 -o main -I/usr/local/include main.cpp -lpapi
generate: generate.cpp
	g++ -std=c++11 -o generate generate.cpp
compare: compare.cpp
	g++ -std=c++11 -o compare compare.cpp
form: form.cpp
	g++ -std=c++11 -o form form.cpp

test: compare main
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0 0
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1 0
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2 0
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0 1
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1 1
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2 1
	./compare "./test/C0.dat" "./test/C.dat"
	./compare "./test/C1.dat" "./test/C.dat"
	./compare "./test/C2.dat" "./test/C.dat"
	rm -f ./test/C0.dat ./test/C1.dat ./test/C2.dat

report: generate main
	for number in 1000 2000 3000 4000 5000; do \
	./generate f $$number $$number A.dat ;\
	./generate f $$number $$number B.dat ;\
	for i in 0 1; do \
	./main A.dat B.dat C.dat 0 $$i;\
	./main A.dat B.dat C.dat 1 $$i;\
	./main A.dat B.dat C.dat 2 $$i;\
	done \
	done

plot: form
	./form cyc.txt ./plot/cycplot.txt
	./form l1.txt ./plot/l1plot.txt
	./form l2.txt ./plot/l2plot.txt
	./form time.txt ./plot/timeplot.txt
	./form flops.txt ./plot/flopsplot.txt
	gnuplot plotset

clean:
	rm -f *.dat *.txt *.o *.svg print compare generate main form
