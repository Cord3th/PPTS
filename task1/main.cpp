#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
using namespace std;

template <typename Atype, typename Btype, typename Ctype>
void mulmatr(fstream &a, fstream &b, fstream &c, int &r){
	uint64_t na, ma, nb, mb, i, j, k;
	a.read((char *) &na, sizeof(na));
	a.read((char *) &ma, sizeof(ma));
	b.read((char *) &nb, sizeof(nb));
	b.read((char *) &mb, sizeof(mb));
	if (ma != nb) {
		cout << "cant multiply matrs because of size\n";
		return;
	}
	Atype **matra = new Atype*[na];
	Btype **matrb = new Btype*[nb];
	Ctype **matrc = new Ctype*[na];
	for (i = 0; i < na; ++i) {
		matra[i] = new Atype[ma];
		matrc[i] = new Ctype[mb];
	}
	for (i = 0; i < nb; ++i) {
		matrb[i] = new Btype[mb];
	}

	for (i = 0 ; i < na; ++i) {
		for(j = 0; j < ma; ++j) {
			a.read((char *) &matra[i][j], sizeof(matra[i][j]));
		}
	}
	for (i = 0 ; i < nb; ++i){
		for (j = 0; j < mb; ++j) {
			b.read((char *) &matrb[i][j], sizeof(matrb[i][j]));
		}
	}
	for( i = 0 ; i < na; ++i) {
		for(j = 0 ; j < mb; ++j) {
			matrc[i][j] = 0;
		}
	}

	ofstream dat;
	dat.open("dat.txt", ios::out | ios::app);

	clock_t time;
	switch (r) {
		case 0 :
			time = clock();
			for(i = 0 ; i < na; ++i) {
				for( j = 0 ; j < mb; ++j) {
					for ( k = 0; k < ma; ++k) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 1 :
			time = clock();
			for(i = 0 ; i < na; ++i) {
				for( k = 0 ; k < ma; ++k) {
					for ( j = 0; j < mb; ++j) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 2 :
			time = clock();
			for(k = 0 ; k < ma; ++k) {
				for( i = 0 ; i < na; ++i) {
					for ( j = 0; j < mb; ++j) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 3 :
			time = clock();
			for(j = 0 ; j < mb; ++j) {
				for( i = 0 ; i < na; ++i) {
					for ( k = 0; k < ma; ++k) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 4 :
			time = clock();
			for(j = 0 ; j < mb; ++j) {
				for( k = 0 ; k < ma; ++k) {
					for ( i = 0; i < na; ++i) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		case 5 :
			time = clock();
			for(k = 0 ; k < ma; ++k) {
				for( j = 0 ; j < mb; ++j) {
					for ( i = 0; i < na; ++i) {
						matrc[i][j] += matra[i][k] * matrb[k][j];
					}
				}
			}
			time -= clock();
			break;
		default : 
			cout << "error: mode" << r << endl;
			return;
	}
	dat << r << ' ' << fixed << setprecision(6) << ((double) -time)/CLOCKS_PER_SEC << endl;
	char type;
	if ( sizeof(Ctype) == sizeof(double)) {
		type = 'd';
	} else {
		type = 'f';
	}

	c.write((char*) &type, sizeof(type));
	c.write((char*) &na, sizeof(na));
	c.write((char*) &mb, sizeof(mb));
	for (i = 0 ; i < na; ++i) {
		for ( j = 0; j < mb; ++j) {
			c.write((char*) &matrc[i][j], sizeof(matrc[i][j]));
		}
	}


	for (i = 0; i < na; ++i) {
		delete []matra[i];
		delete []matrc[i];
	}
	for (i = 0; i < nb; ++i) {
		delete []matrb[i];
	}
	delete []matra;
	delete []matrb;
	delete []matrc;
}

int main(int argc, char **argv) {
	//ijk ikj kij jik jki kji
	
	if (argc != 5) {
		cout << "format: A.dat B.dat input.dat 0/1/2/3/4/5" << endl;
		return 0;
	}

	fstream a,b,c;
	a.open(argv[1], ios::in | ios::binary);
	b.open(argv[2], ios::in | ios::binary);
	char ta, tb;

	c.open(argv[3], ios::out | ios::binary);
	int r;
	sscanf(argv[4], "%d", &r);
	a.read(&ta, sizeof(ta));
	b.read(&tb, sizeof(tb));
	if (ta == 'f' && tb == 'f') {
		mulmatr <float, float, float>(a, b, c, r);
	} else if(ta == 'f' && tb == 'd'){
		mulmatr <float, double, double>(a, b, c, r);
	} else if(ta == 'd' && tb == 'f') {
		mulmatr <double, float, double>(a, b, c, r);
	} else if ( ta == 'd' && tb == 'd') {
		mulmatr <double, double, double>(a, b, c, r);
	} else {
		cout << "error with types of matrs" << endl;
	}
	a.close();
	b.close();
	c.close();

	return 0;
}