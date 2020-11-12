#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char **argv) {
    int rank, size;
    const int data_tag = 1, time_tag = 2;
    long long first, last, temp;
    double time_start, time_finish;
    MPI_Status status;
    vector<double> time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sscanf(argv[1], "%llu", &first);
    sscanf(argv[2], "%llu", &last);

    time_start = MPI_Wtime();
    vector<bool> sqrt_primes((int) sqrt(last + 1), true);

    sqrt_primes[0] = sqrt_primes[1] = false;
    for (temp = 2; temp * temp <= last; temp++) {
        if (sqrt_primes[temp]) {
            for (long long j = temp * temp; j * j <= last; j += temp) {
                sqrt_primes[j] = false;
            }
        }
    }

    time_finish = MPI_Wtime();
    time.push_back(time_finish - time_start);

    long long i_first = (long long) (rank - 1)
                       * (last - max(first, temp) + 1)
                       / (size - 1)
                       + max(first, temp);
    long long i_last = (long long) rank
                     * (last - max(first, temp) + 1)
                     / (size - 1)
                     + max(first, temp) - 1;

    if (rank) {
        time_start = MPI_Wtime();
        vector<bool> i_primes(i_last - i_first + 1, true);
        long long i;

        for (long long j = 2; j * j <= last; j++) {
            for (i = (i_first / j + 1 * (i_first % j != 0)) * j; i <= min(i_last, last); i += j) {
                i_primes[i - i_first] = false;
            }
        }

        for (i = i_first; i <= min(i_last, last); i++) {
            if (i_primes[i - i_first]) {
                MPI_Send(&i, 1, MPI_LONG_LONG, 0, data_tag, MPI_COMM_WORLD);
            }
        }

        i = -1;
        MPI_Send(&i, 1, MPI_LONG_LONG, 0, data_tag, MPI_COMM_WORLD);

        time_finish = MPI_Wtime();
        double i_time = time_finish - time_start;
        MPI_Send(&i_time, 1, MPI_DOUBLE, 0, time_tag, MPI_COMM_WORLD);
    } else {
        int k = 0;
        double sum_time = 0, max_time = 0;
        long long prime_count = 0, size_count = 0;
        vector<long long> remain_primes;

        while (size_count < (size - 1)) {
            long long tmp;

            MPI_Recv(&tmp, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, data_tag, MPI_COMM_WORLD, &status);
            
            if (tmp != -1) {
                remain_primes.push_back(tmp);
                prime_count++;
            } else {
                double tmptime;
                MPI_Recv(&tmptime, 1, MPI_DOUBLE, MPI_ANY_SOURCE, time_tag, MPI_COMM_WORLD, &status);
                time.push_back(tmptime);
                size_count++;
            }
        }

        ofstream out;
        out.open(argv[3]);

        for (long long i = 0; i < prime_count; i++) {
            out << remain_primes[i] << '\n';
        }

        for (temp = first; temp * temp <= last; temp++) {
            if (sqrt_primes[temp]) {
                out << temp << '\n';
                prime_count++;
            }
        }

        cout << "There are " << prime_count << " primes" << endl;


        for (k = 1; k < size; k++) {
            sum_time += time[k - 1];
            if (time[k - 1] > max_time) {
                max_time = time[k - 1];
            }
        }

        cout << "Overall time: " << sum_time << endl
             << "Maximal single process time: " << max_time << endl;
    }

    MPI_Finalize();

    return 0;
}
