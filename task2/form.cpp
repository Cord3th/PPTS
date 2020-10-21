#include <fstream>
#include <iostream>
using namespace std;


int main(int argc, char **argv) {
	if (argc != 3) {
		cerr << "Format: input output" << endl;
		return 0;
	}
	ifstream a;
	ofstream b;
	b.open(argv[2]);
	a.open(argv[1]);
	long long n;
	double ar[5][3];
	while (a >> n) {
		long long r;
		a >> r;
		double val;
		a >> val;
		ar[n / 1000 - 1][r] = val;
	}
	for (int i = 0; i < 5; i++) {
		b << (i + 1) * 1000 << ' ' << ar[i][0] << ' ' << ar[i][1] << ' ' << ar[i][2] << endl;
	}
	return 0;
}
