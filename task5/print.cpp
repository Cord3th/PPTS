#include <fstream>
#include <iostream>
#include <iomanip>


using namespace std;
int main(int argc, char** argv) {

	if (argc != 3) {
		cout << "Format: input.dat otput.txt" << endl;
		return 1;
	}

	const int w = 8;
	fstream input_file;
	ofstream output_file(argv[2]);
	char type;
	uint64_t n, m;

	input_file.open(argv[1], ios::binary | ios::in);
	input_file.read(&type, sizeof(type));
	input_file.read((char *) &n, sizeof(n));
	input_file.read((char *) &m, sizeof(m));
	for (uint64_t i = 0; i < n; i++) {
		for (uint64_t j = 0; j < m; j++) {
			if (type == 'f') {
				float tmp;
				input_file.read((char *) &tmp, sizeof(tmp));
				output_file << setw(w) << tmp << ' ';
			} else if (type == 'd') {
				double tmp;
				input_file.read((char *) &tmp, sizeof(tmp));
				output_file << setw(w) << tmp << ' ';
			} else {
				cerr << "Matrix type error" << endl;
				return 2;
			}
		}
		output_file << endl;
	}

	input_file.close();
	output_file.close();

	return 0;
}
