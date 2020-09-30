#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char ** argv) {
	if ( argc != 5) {
		cout << "format: Type(f/d) n m  fout.dat" <<endl;
		return 2;
	}

	srand(time(0));
	char *type = argv[1];
	uint64_t n,m;
	sscanf(argv[2], "%lu", &n);
	sscanf(argv[3], "%lu", &m);
	string filename(argv[4]);
	fstream file;
	file.open(filename, ios::out | ios::binary | ios::trunc);
	file.write(type, sizeof(char));
	file.write((char *) &n, sizeof(n));
	file.write((char *) &m, sizeof(m));
	uint64_t i, j;
	for ( i = 0 ; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			if ( *type == 'f') {
				float temp = ((float) rand())/rand() - ((float) rand())/rand();
				file.write((char *) &temp, sizeof(temp));
			} else if ( *type == 'd') {
				double temp = ((double) rand())/rand() - ((double) rand())/rand();
				file.write((char *) &temp, sizeof(temp));
			} else {
				cout << "format: Type(f/d) n m  fout.dat" << endl;
				return 1;
			}
		}
	}
	file.close();



	return 0;
}