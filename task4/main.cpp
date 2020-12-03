#include "mpi.h"
#include <iostream>
#include <stdint.h>
#include <fstream>

using namespace std;


int main(int argc, char **argv) {
	if (argc != 4) {
		cerr << "Format: A.dat b.dat c.dat" << endl;
		return 0;
	}

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	const int data_tag = 1, b_tag = 2, c_tag = 3;
	double time_start, time_finish, time, max_time, sum_time = 0;
	fstream b_input_file(argv[2]);
	uint64_t b_m, b_n, n, m;
	char type;

	b_input_file.read(&type, sizeof(type));
	b_input_file.read((char *) &b_m, sizeof(b_m));
	b_input_file.read((char *) &b_n, sizeof(b_n));

	if (rank == 0) {
		fstream input_file(argv[1]);
		input_file.read(&type, sizeof(type));
		input_file.read((char *) &m, sizeof(m));
		input_file.read((char *) &n, sizeof(n));

		cout << "m = " << m <<" n = " << n << " size = " << size << endl;

		int sizes[2];
		sizes[1] = n;
		sizes[0] = m;

		if (m >= n) {
			long long begin = (long long) 0;
			long long end = (long long) m / size - 1;

			double a[end + 1][n];

			for (int i = begin; i <= end; i++) {
				for (int j = 0; j < n; j++) {
					input_file.read((char *) &a[i][j], sizeof(a[i][j]));
				}
			}

    		for (int k = 1; k < size; k++) { // k = rank
    			MPI_Send(&sizes[0], 2, MPI_INT, k, data_tag, MPI_COMM_WORLD);

    			begin = (long long) k * m / size ;
    			end = (long long) (k + 1) * m / size - 1;
    			double temp[end - begin + 1][n];

    			for (int i = begin; i <= end; i++) {
    				for (int j = 0; j < n; j++) {
    					input_file.read((char *) &temp[i - begin][j],
										sizeof(temp[i - begin][j]));
    				}
    				MPI_Send(&temp[i-begin][0], n, MPI_DOUBLE, k, data_tag, MPI_COMM_WORLD);
    			}
    		}

    		double b[b_m], c[m];

    		for (int i = 0; i < b_m; i++) {
    			b_input_file.read((char *) &b[i], sizeof(b[i]));
    		}

    		time_start = MPI_Wtime();

    		begin = (long long) 0;
    		end = (long long) m / size - 1;

    		for (int i = begin; i <= end; i++) {
    			c[i - begin] = 0;
    			for (int j = 0; j < n; j++) {
    				c[i - begin] += a[i - begin][j] * b[j];
    			}
    		}

    		time_finish = MPI_Wtime();
    		time = time_finish - time_start;

    		for (int i = 1; i < size; i++) {
    			begin = (long long) i * m / size ;
    			end = (long long) (i + 1) * m / size - 1;
    			MPI_Recv(&c[begin], end - begin + 1, MPI_DOUBLE, i, c_tag, MPI_COMM_WORLD, &status);
    		}

    		fstream c_output_file;
    		c_output_file.open(argv[3], ios::out | ios::binary);
    		type = 'd';
    		c_output_file.write(&type, sizeof(type));
    		c_output_file.write((char *) &m, sizeof(m));
    		n = 1;
    		c_output_file.write((char *) &n, sizeof(n));
    		for (int i = 0; i < m; i++) {
    			c_output_file.write((char *) &c[i], sizeof(c[i]));
    		}
    	} else {
    		double **a, b[b_m];
			a = new double *[n];

    		for (int i = 0; i < m; i++) {
				a[i] = new double[n];
    			for (int j = 0; j < n; j++) {
    				input_file.read((char *) &a[i][j], sizeof(a[i][j]));
    			}
    		}

    		for (int i = 0; i < b_m; i++) {
    			b_input_file.read((char *) &b[i], sizeof(b[i]));
    		}

    		for (int k = 1; k < size; k++) {
    			MPI_Send(&sizes[0], 2, MPI_INT, k, data_tag, MPI_COMM_WORLD);

    			long long begin = (long long) k * n / size;
    			long long end = (long long) (k + 1) * n / size - 1;
    			double temp[end - begin + 1][m];

    			for (int j = begin; j <= end; j++) {
    				for (int i = 0; i < m; i++) {
    					temp[j - begin][i] = a[i][j];
    				}
    				MPI_Send(&temp[j - begin][0], m, MPI_DOUBLE, k, data_tag, MPI_COMM_WORLD);
    			}

    			MPI_Send(&b[begin], end - begin + 1, MPI_DOUBLE, k, b_tag, MPI_COMM_WORLD);
    		}

    		double c[m], c1[m];
    		time_start = MPI_Wtime();
    		long long begin = (long long) 0;
    		long long end = (long long) n / size - 1;

    		for (int i = 0; i < m; i++) {
    			c[i] = 0;
    			for (int j = begin; j <= end; j++) {
    				c[i] += a[i][j - begin] * b[j - begin];
    			}
    		}

    		time_finish = MPI_Wtime();
    		sum_time = time_finish - time_start;
    		time = sum_time;

    		MPI_Reduce(c, c1, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    		fstream c_output_file;

	    	c_output_file.open(argv[3], ios::out | ios::binary);
	    	type = 'd';
	    	c_output_file.write(&type, sizeof(type));
	    	c_output_file.write((char *) &m, sizeof(m));
	    	n = 1;
	    	c_output_file.write((char *) &n, sizeof(n));

	    	for (int i = 0; i < m; i++) {
	    		c_output_file.write((char *) &c1[i], sizeof(c[i]));
			}
    	}

    } else {
    	int sizes[2];

    	MPI_Recv(&sizes[0], 2, MPI_INT, 0, data_tag, MPI_COMM_WORLD, &status);

    	m = sizes[0];
    	n = sizes[1];

    	if (m >= n) {
			long long begin = (long long) rank * m / size ;
			long long end = (long long) (rank + 1) * m / size - 1;
    		double b[b_m], a[end - begin + 1][n], c[end - begin + 1];

    		for (int i = 0; i < b_m; i++) {
    			b_input_file.read((char *) &b[i], sizeof(b[i]));
    		}


    		for (int i = begin; i <= end; i++) {
    			MPI_Recv(&a[i - begin][0], n, MPI_DOUBLE, 0, data_tag, MPI_COMM_WORLD, &status);
    		}

    		time_start = MPI_Wtime();

    		for (int i = begin; i <= end; i++) {
    			c[i - begin] = 0;
    			for (int j = 0; j < n; j++) {
    				c[i - begin] += a[i - begin][j] * b[j];
    			}
    		}

    		time_finish = MPI_Wtime();

    		MPI_Send(&c[0], end - begin  + 1, MPI_DOUBLE, 0, c_tag, MPI_COMM_WORLD);

    		time = time_finish - time_start;
    	} else {
    		long long begin = (long long) rank * n / size ;
    		long long end = (long long) (rank + 1) * n / size - 1;
    		double a[end - begin + 1][m], b[end - begin + 1], c[m];

    		for (int i = begin; i <= end; i++) {
    			MPI_Recv(&a[i - begin][0], m, MPI_DOUBLE, 0, data_tag, MPI_COMM_WORLD, &status);
    		}

    		MPI_Recv(&b[0], end - begin + 1, MPI_DOUBLE, 0, b_tag, MPI_COMM_WORLD, &status);

			time_start = MPI_Wtime();

    		for (int i = 0; i < m; i++) {
    			c[i] = 0;
    			for (int j = begin; j <= end; j++) {
    				c[i] += b[j - begin] * a[j - begin][i];
    			}
    		}

    		time_finish = MPI_Wtime();

    		time = time_finish - time_start;
    		MPI_Reduce(c, 0, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	}
    }

    MPI_Reduce(&time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0 ) {
    	cout << "Overall time :" << sum_time << endl;
    	cout << "Maximal single process time :" << max_time << endl;
    }

    MPI_Finalize();
    return 0;
}
