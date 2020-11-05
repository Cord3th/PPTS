#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double time_start, time_finish;
    long long  begin, end;
    int tagdata = 1;
    int tagtime = 2;
    sscanf(argv[1], "%llu", &begin);
    sscanf(argv[2], "%llu", &end);
    vector<bool> prime((int) sqrt(end + 1), true);

    MPI_Status status;
    prime[0] = prime[1] = false;
    long long  temp;
    for (temp = 2; temp * temp <= end; temp++) {
        if (prime[temp]) {
            for (long long j = temp * temp; j * j <= end; j += temp) {
                prime[j] = false;
            }
        }

    }

    long long  ibegin = (long long) (rank - 1) * (end - max(begin,temp) + 1)
                        / (size - 1) + max(begin,temp);
    long long  iend = (long long ) rank * (end - max(begin,temp) + 1)
                      / (size - 1) + max(begin,temp) - 1;

    if (rank != 0) {
        time_start = MPI_Wtime();
        vector<bool> prime1(iend - ibegin + 1, true);
        for (long long j = 2; j * j <= end; j++) {
            if (prime[j]) {
                for (long long i = (ibegin / j + 1 * (ibegin % j != 0)) * j;
                     i <= min(iend, end); i += j) {
                    prime1[i - ibegin] = false;
                }
            }

        }
        for (long long i = ibegin; i <= min(iend, end); i++) {
            if (prime1[i - ibegin]) {
                MPI_Send(&i, 1 , MPI_LONG_LONG, 0, tagdata, MPI_COMM_WORLD);
            }
        }
        long long int i = -1;
        MPI_Send(&i, 1, MPI_LONG_LONG, 0, tagdata, MPI_COMM_WORLD);
        time_finish = MPI_Wtime();
        double time = time_finish - time_start;
        MPI_Send(&time, 1, MPI_DOUBLE, 0, tagtime, MPI_COMM_WORLD);
    } else {
        std::vector<long long> simple;
        long long count = 0;
        long long countfin = 0;
        std::vector<double> time;
        int k = 0;
        while(1) {
            long long  tmp ;

            MPI_Recv(&tmp, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, tagdata, MPI_COMM_WORLD, &status);
            if (tmp != -1) {
                simple.push_back(tmp);
                count++;
            } else {
                double tmptime;
                MPI_Recv(&tmptime, 1, MPI_DOUBLE ,MPI_ANY_SOURCE, tagtime, MPI_COMM_WORLD, &status);
                time.push_back(tmptime);
                countfin++;
                if (countfin == (size - 1)) {
                    break;
                }
            }
        }

        ofstream out;
        out.open(argv[3]);
        for (temp = begin; temp * temp <= end; temp++) {
            if (prime[temp]) {
                out << temp << ' ';
                count++;
            }

        }

        cout << count << endl;
        for (long long i = 0; i < count; i++) {
            out << simple[i] << '\n';
        }
        double sumtime = 0;
        double maxtime = 0;
        for (k = 1; k < size; k++) {
            sumtime += time[k-1];
            if (time[k-1] > maxtime) {
                maxtime = time[k-1];
            }
        }
        cout << "all time: " << sumtime << endl << "max time: " << maxtime << endl;
    }
    MPI_Finalize();
    return 0;
}
