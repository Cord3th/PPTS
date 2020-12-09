#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv) {
    int rank, size, first = atoi(argv[1]), last = atoi(argv[2]),
        sqrt_last = sqrt(last), i_min = max(first - 1, sqrt_last);
    const int data_tag = 1, time_tag = 2;
    double time_start, time_finish;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Barrier(MPI_COMM_WORLD);

    vector<int> sqrt_primes(sqrt_last);
    time_start = MPI_Wtime();

    sqrt_primes[0] = 0;
    for (size_t i = 1; i < sqrt_last; i++) {
        sqrt_primes[i] = i + 1;
    }

    for (size_t i = 1; i < sqrt_last; i++) {
        if (sqrt_primes[i]) {
            for (int j = i + sqrt_primes[i]; j < sqrt_last; j += sqrt_primes[i]) {
                sqrt_primes[j] = 0;
            }
        }
    }

    int step = (last - i_min) / size
               + ((last - i_min) % size != 0),
        i_first = i_min + rank * step + 1;

    int *i_primes = new int[step + 1],
        i_max = min(i_first + step, last);
    for (int i = 0; i < step; i++) {
        if ((i + i_first) < i_max) {
            i_primes[i] = i + i_first;
        } else {
            i_primes[i] = 0;
        }
    }

    for (int i = 0; i < sqrt_last; i++) {
        if (sqrt_primes[i]) {
            for (int j = (i_first / sqrt_primes[i] + 1 * (i_first % sqrt_primes[i] != 0)) * sqrt_primes[i] - 1;
                j < i_max; j += sqrt_primes[i]) {
                    i_primes[j - i_first + 1] = 0;
                }
        }
    }

    time_finish = MPI_Wtime();
    if (rank) {
        MPI_Send(i_primes, step, MPI_INT, 0, data_tag, MPI_COMM_WORLD);
        double time = time_finish - time_start;
        MPI_Send(&time, 1, MPI_DOUBLE, 0, time_tag, MPI_COMM_WORLD);
    } else {
        //ofstream out;
        //out.open(argv[3]);
        long long prime_count = 0;
        double sum_time = time_finish - time_start, max_time = sum_time;
        for (int i = first - 1; i < sqrt_last; i++) {
            if (sqrt_primes[i]) {
                //out << sqrt_primes[i] << '\n';
                prime_count++;
            }
        }
        for (int i = 0; i < step; i++) {
            if (i_primes[i]) {
                //out << i_primes[i] << '\n';
                prime_count++;
            }
        }
        for (int i = 0; i < size - 1; i++) {
			MPI_Recv(i_primes, step, MPI_INT, MPI_ANY_SOURCE, data_tag, MPI_COMM_WORLD, &status);
			for (int j = 0; j < step; j++) {
                if (i_primes[j]) {
                    prime_count++;
                    //out << i_primes[j] << '\n';
                }
            }
		}
        for (int i = 0; i < size - 1; i++) {
			double time;
			MPI_Recv(&time, 1, MPI_DOUBLE, MPI_ANY_SOURCE, time_tag, MPI_COMM_WORLD, &status);
			if (time > max_time)
				max_time = time;
			sum_time += time;
		}
        cout << prime_count << endl;
        /*cout << "There are " << prime_count << " primes" << endl
             << "Overall time: " << sum_time << endl
             << "Maximal single process time: " << max_time << endl;*/
    }
    delete[] i_primes;
    MPI_Finalize();

    return 0;
}
