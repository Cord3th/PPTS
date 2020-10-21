#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <papi.h>

#define NUM_EVENT 3
#define THRESHOLD 100000
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }

using namespace std;

void mulmatr(fstream &a, fstream &b, fstream &c, int mode, int m) {
	int retval, EventSet = PAPI_NULL,
	 	event_codes[NUM_EVENT] = {PAPI_TOT_CYC, PAPI_L1_DCM, PAPI_L2_DCM};
    char errstring[PAPI_MAX_STR_LEN];
	float ptime, rtime, mflops;
    long long flpops, values[NUM_EVENT];

	uint64_t A_rows, A_cols, B_rows, B_cols,
	 		 i, j, k, n, i1, j1, k1, blocksize = 32,
			 myblock = 72;

	ofstream cyc, l1, l2, time, flops;
	cyc.open("cyc.txt", ios::out | ios::app);
	l1.open("l1.txt", ios::out | ios::app);
	l2.open("l2.txt", ios::out | ios::app);
	time.open("time.txt", ios::out | ios::app);
	flops.open("flops.txt", ios::out | ios::app);

	a.read((char *) &A_rows, sizeof(A_rows));
	a.read((char *) &A_cols, sizeof(A_cols));
	b.read((char *) &B_rows, sizeof(B_rows));
	b.read((char *) &B_cols, sizeof(B_cols));

	if (A_cols != B_rows) {
		cerr << "Size incompatibility\n";
		return;
	}

	n = A_rows;

	float **A = new float*[A_rows];
	float **B = new float*[B_rows];
	float **C = new float*[A_rows];

	for (i = 0; i < A_rows; i++) {
		A[i] = new float[A_cols];
		C[i] = new float[B_cols];
	}

	for (i = 0; i < B_rows; i++) {
		B[i] = new float[B_cols];
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

    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
		fprintf(stderr, "Error: %s\n", errstring);
    	exit(1);
    }

	switch (mode) {
		case 0:
			if (m == 0) {
				if ((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_start(EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}

			for (i = 0; i < A_rows; i += blocksize) {
				for (j = 0; j < B_cols; j += blocksize) {
					for (k = 0; k < A_cols; k += blocksize) {
						//block mul
						uint64_t imin = min(n, i + blocksize),
								 jmin = min(n, j + blocksize),
								 kmin = min(n, k + blocksize);
						for (i1 = i; i1 < imin; i1++) {
							for (j1 = j; j1 < jmin; j1++) {
								for (k1 = k; k1 < kmin; k1++) {
									C[i1][j1] += A[i1][k1] * B[k1][j1];
								}
							}
						}
					}
				}
			}
			if (m == 0) {
			    if ((retval = PAPI_stop(EventSet,values)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
				if ((retval = PAPI_remove_events(EventSet,event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

			    if ((retval = PAPI_destroy_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}
			break;

		case 1:
			if (m == 0) {
				if ((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_start(EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}
			for (i = 0; i < A_rows; i += blocksize) {
				for (j = 0; j < B_cols; j += blocksize) {
					for (k = 0; k < A_cols; k += blocksize) {
						//block mul
						uint64_t imin = min(n, i + blocksize),
								 jmin = min(n, j + blocksize),
								 kmin = min(n, k + blocksize);
						for (i1 = i; i1 < imin; i1++) {
							for (k1 = k; k1 < kmin; k1++) {
								for (j1 = j; j1 < jmin; j1++) {
									C[i1][j1] += A[i1][k1] * B[k1][j1];
								}
							}
						}
					}
				}
			}

			if (m == 0) {
			    if ((retval = PAPI_stop(EventSet,values)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
				if ((retval = PAPI_remove_events(EventSet,event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

			    if ((retval = PAPI_destroy_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}
			break;

		case 2:
			blocksize = myblock;
			if (m == 0) {
				if ((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

				if ((retval = PAPI_start(EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}
			for (i = 0; i < A_rows; i += blocksize) {
				for (j = 0; j < B_cols; j += blocksize) {
					for (k = 0; k < A_cols; k += blocksize) {
							//block mul
						uint64_t imin = min(n, i + blocksize),
								 jmin = min(n, j + blocksize),
								 kmin = min(n, k + blocksize);
						for (i1 = i; i1 < imin; i1++) {
							for (k1 = k; k1 < kmin; k1++) {
								for (j1 = j; j1 < jmin; j1++) {
									C[i1][j1] += A[i1][k1] * B[k1][j1];
								}
							}
						}
					}
				}
			}
			if (m == 0) {
			    if ((retval = PAPI_stop(EventSet,values)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
				if ((retval = PAPI_remove_events(EventSet,event_codes, NUM_EVENT)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}

			    if ((retval = PAPI_destroy_eventset(&EventSet)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			} else {
				if ((retval = PAPI_flops(&rtime, &ptime, &flpops, &mflops)) != PAPI_OK) {
					ERROR_RETURN(retval);
				}
			}
			break;
		default :
			cerr << "Error: mode" << mode << endl;
			return;
	}
	if (m == 0) {
		cyc << n << ' ' << mode << ' ' << values[0] << endl;
		l1 << n << ' ' << mode  << ' ' << values[1] << endl;
		l2 << n << ' ' << mode  << ' ' << values[2] << endl;
	} else {
		time << n << ' ' << mode << ' ' << ptime << endl;
		flops << n << ' ' << mode << ' ' << flpops << endl;
	}

	cyc.close();
	l1.close();
	l2.close();
	time.close();
	flops.close();

    PAPI_shutdown();

	char type ='f';

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
	//ijk == 0 , ikj == 1

	if (argc != 6) {
		cerr << "Format: A.dat B.dat input.dat  mode (0/1/2) m(0/1)" << endl;
		return 0;
	}

	fstream a,b,c;
	a.open(argv[1], ios::in | ios::binary);
	b.open(argv[2], ios::in | ios::binary);
	char ta, tb;

	c.open(argv[3], ios::out | ios::binary);
	int mode;
	sscanf(argv[4], "%d", &mode);
	int m;
	sscanf(argv[5], "%d", &m);
	a.read(&ta, sizeof(ta));
	b.read(&tb, sizeof(tb));
	mulmatr(a, b, c, mode, m);
	a.close();
	b.close();
	c.close();

	return 0;
}
