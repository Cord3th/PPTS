#include <fstream>
#include <iostream>
#include <iomanip>
const int w = 8;

using namespace std;
int main(int argc, char** argv) {
	if (argc != 3) {
		cout << "format: iput.dat otput.txt" << endl;
		return 1;
	}
	fstream fin;
	fin.open(argv[1], ios::binary | ios::in);
	ofstream fout(argv[2]);
	char type;
	fin.read(&type, sizeof(type));
	uint64_t n, m, i, j;
	fin.read((char *) &n, sizeof(n));
	fin.read((char *) &m, sizeof(m));
	for(i = 0 ; i < n; ++i) {
		for (j = 0 ; j < m; ++j) {
			if (type == 'f') {
				float tmp;
				fin.read((char *) &tmp, sizeof(tmp));
				fout << setw(w) << tmp << ' ';
			} else if (type == 'd'){
				double tmp;
				fin.read((char *) &tmp, sizeof(tmp));
				fout << setw(w) << tmp << ' ';
			} else {
				cout << "error matr" << endl;
				return 2;
			}
		}
		fout << endl;
	}
	fin.close();
	fout.close();


	return 0;
}