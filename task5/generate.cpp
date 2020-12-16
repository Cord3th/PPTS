#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char ** argv) {
	if (argc != 5) {
		cerr << "Format: Type (f/d) n m fout.dat" <<endl;
		return 2;
	}

	char *type = argv[1];
	uint64_t n, m;
	string filename(argv[4]);
	fstream file;

	srand(time(0));

	sscanf(argv[2], "%lu", &n);
	sscanf(argv[3], "%lu", &m);

	file.open(filename, ios::out | ios::binary | ios::trunc);
	file.write(type, sizeof(char));
	file.write((char *) &n, sizeof(n));
	file.write((char *) &m, sizeof(m));

	for (uint64_t i = 0 ; i < n; i++) {
		for (uint64_t j = 0; j < m; j++) {
			if (*type == 'f') {
				float temp = ((float) rand() / (RAND_MAX));
				file.write((char *) &temp, sizeof(temp));
			} else if (*type == 'd') {
				double temp = ((double) rand() / (RAND_MAX));
				file.write((char *) &temp, sizeof(temp));
			} else {
				cerr << "Format: Type (f/d) n m fout.dat" << endl;
				return 1;
			}
		}
	}
	file.close();

	return 0;
}
