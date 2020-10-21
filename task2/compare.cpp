#include <iostream>
#include <fstream>
#include <cfloat>
#include <cmath>

using namespace std;

const double eps = DBL_EPSILON;


int main(int argc, char ** argv) {
	if (argc != 3) {
		cerr << "Format : A.dat B.dat" << endl;
		return 1;
	}
	fstream f1,f2;
	char type1, type2;
	uint64_t n, n1, m1, m, i, j;
	f1.open(argv[1], ios::in | ios::binary);
	f2.open(argv[2], ios::in | ios::binary);
	f1.read((char *) &type1, sizeof(type1));
	f2.read((char *) &type2, sizeof(type2));
	if (type1 != type2) {
		cerr << "Type diff\n";
		f1.close();
		f2.close();
		return 0;
	}
	f1.read((char *) &n, sizeof(n));
	f2.read((char *) &n1, sizeof(n1));
	if (n != n1) {
		cerr << "Size diff\n";
		f1.close();
		f2.close();
		return 0;
	}
	f1.read((char *) &m, sizeof(m));
	f2.read((char *) &m1, sizeof(m1));
	if (m != m1) {
		cerr << "Size diff\n";
		f1.close();
		f2.close();
		return 0;
	}
	if (type1 == 'f') {
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				float tmp1, tmp2;
				f1.read((char *) &tmp1, sizeof(tmp1));
				f2.read((char *) &tmp2, sizeof(tmp2));
				if (abs(tmp1 - tmp2) > eps) {
					cerr << "Elements diff\n";
					f1.close();
					f2.close();
					return 0;
				}
			}
		}
	} else {
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				double tmp1, tmp2;
				f1.read((char *) &tmp1, sizeof(tmp1));
				f2.read((char *) &tmp2, sizeof(tmp2));
				if (abs(tmp1 - tmp2) > eps) {
					cerr << "Elements diff\n";
					f1.close();
					f2.close();
					return 0;
				}
			}
		}
	}
	cout << "The matrices are equal\n";
	return 0;

}
