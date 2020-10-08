#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <iomanip>

using namespace std;

template <typename Mtype>
void mulmatr(fstream &a, fstream &b, fstream &c, int &mode){
	uint64_t A_rows, A_cols, B_rows, B_cols, i, j, k;
	a.read((char *) &A_rows, sizeof(A_rows));
	a.read((char *) &A_cols, sizeof(A_cols));
	b.read((char *) &B_rows, sizeof(B_rows));
	b.read((char *) &B_cols, sizeof(B_cols));
	if (A_cols != B_rows) {
		cerr << "Size incompatibility \n";
		return;
	}
	Mtype **A = new Mtype*[A_rows];
	Mtype **B = new Mtype*[B_rows];
	Mtype **C = new Mtype*[A_rows];
	for (i = 0; i < A_rows; i++) {
		A[i] = new Mtype[A_cols];
		C[i] = new Mtype[B_cols];
	}
	for (i = 0; i < B_rows; i++) {
		B[i] = new Mtype[B_cols];
	}

	for (i = 0; i < A_rows; i++) {
		for (j = 0; j < A_cols; j++) {
			a.read((char *) &A[i][j], sizeof(A[i][j]));
		}
	}
	for (i = 0; i < B_rows; i++){
		for (j = 0; j < B_cols; j++) {
			b.read((char *) &B[i][j], sizeof(B[i][j]));
		}
	}
	for (i = 0; i < A_rows; i++) {
		for (j = 0; j < B_cols; j++) {
			C[i][j] = 0;
		}
	}

	ofstream dat;
	dat.open("dat.txt", ios::out | ios::app);

	clock_t time;
	switch (mode) {
		case 0:
			time = clock();
			for (i = 0; i < A_rows; i++) {
				for (j = 0; j < B_cols; j++) {
					for (k = 0; k < A_cols; k++) {
						C[i][j] += A[k][i] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 1:
			time = clock();
			for (i = 0; i < A_rows; i++) {
				for (k = 0; k < A_cols; k++) {
					for (j = 0; j < B_cols; j++) {
						C[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 2:
			time = clock();
			for (k = 0; k < A_cols; k++) {
				for (i = 0; i < A_rows; i++) {
					for (j = 0; j < B_cols; j++) {
						C[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 3:
			time = clock();
			for (j = 0; j < B_cols; j++) {
				for (i = 0; i < A_rows; i++) {
					for (k = 0; k < A_cols; k++) {
						C[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 4:
			time = clock();
			for (j = 0; j < B_cols; j++) {
				for (k = 0; k < A_cols; k++) {
					for (i = 0; i < A_rows; i++) {
						C[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 5:
			time = clock();
			for (k = 0; k < A_cols; k++) {
				for (j = 0; j < B_cols; j++) {
					for (i = 0; i < A_rows; i++) {
						C[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			time -= clock();
			break;
		default :
			cerr << "Error: mode" << mode << endl;
			return;
	}
	dat << mode << ' ' << fixed << setprecision(6) << ((double) -time) / CLOCKS_PER_SEC << endl;
	char type;
	if (sizeof(Mtype) == sizeof(double)) {
		type = 'd';
	} else {
		type = 'f';
	}

	c.write((char*) &type, sizeof(type));
	c.write((char*) &A_rows, sizeof(A_rows));
	c.write((char*) &B_cols, sizeof(B_cols));
	for (i = 0; i < A_rows; i++) {
		for (j = 0; j < B_cols; j++) {
			c.write((char*) &C[i][j], sizeof(C[i][j]));
		}
	}


	for (i = 0; i < A_rows; i++) {
		delete []A[i];
		delete []C[i];
	}
	for (i = 0; i < B_rows; i++) {
		delete []B[i];
	}
	delete []A;
	delete []B;
	delete []C;
}

int main(int argc, char **argv) {
	//ijk ikj kij jik jki kji

	if (argc != 5) {
		cout << "Format: A.dat B.dat input.dat 0/1/2/3/4/5" << endl;
		return 0;
	}

	fstream a, b, c;
	a.open(argv[1], ios::in | ios::binary);
	b.open(argv[2], ios::in | ios::binary);
	char ta, tb;

	c.open(argv[3], ios::out | ios::binary);
	int mode;
	sscanf(argv[4], "%d", &mode);
	a.read(&ta, sizeof(ta));
	b.read(&tb, sizeof(tb));
	if (ta == 'f' && tb == 'f') {
		mulmatr<float>(a, b, c, mode);
	} else if (ta = 'd' && tb == 'd') {
		mulmatr<double>(a, b, c, mode);
	} else {
		cerr << "Matrices type error" << endl;
	}
	a.close();
	b.close();
	c.close();

	return 0;
}
