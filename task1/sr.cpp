#include <fstream>
#include <iostream>

using namespace std;

int main(void) {
	ifstream fin;
	double val[6];
	int count[6];
	int i;
	for(i = 0; i < 6; ++i) {
		val[i] = 0;
		count[i] = 0;
	}
	fin.open("dat.txt");
	int ct;
	while (fin >> ct) {
		double tmp;
		fin >> tmp;
		val[ct] += tmp;
		++count[ct];
	}
	for (i = 0 ; i < 6; ++i) {
		if(count[i] == 0) {
			cout << "err";
			return 1;
		}
		val[i] /= count[i];

	}
	fin.close();

	ofstream fout;
	fout.open("dat.txt");
	for(i = 0 ; i < 6; ++i) {
		fout << i << ' ' << val[i] << endl;
	}
	fout.close();

	return 0;
}