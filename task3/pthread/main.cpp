#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <time.h>


using namespace std;

pthread_mutex_t mutex;

struct smpl {
    long long first;
    long long last;
    int size;
    int rank;
    int prime_count;
    double time;
    std::vector<bool> *primes;
    char *file;
    smpl(int s = 1, int r = 1, int c = 0, long long fr = 0,
         long long l = 1, double t = 0, vector<bool> *p = NULL,
         char *f = 0) {
        size = s;
        rank = r;
        prime_count = c;
        first = fr;
        last = l;
        time = t;
        primes = p;
        file = f;
    }
};

void *worker(void *atr) {
    smpl *a = (smpl*) atr;
    long long i_first = (long long) (a->rank - 1)
                        * (a->last - a->first + 1)
                        / (a->size - 1)
                        + a->first,
    i_last = (long long ) a->rank
             * (a->last - a->first + 1)
             / (a->size - 1)
             + a->first- 1;
    vector<bool> i_primes(i_last - i_first + 1, true);
    timespec time_start, time_finish;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_start);

    for (long long j = 2; j * j <= a->last; j++) {
        if ((*(*a).primes)[j]) {
            for (long long i = (i_first / j + 1 * (i_first % j != 0)) * j;
                i <= min(i_last, a->last); i += j) {
                i_primes[i - i_first] = false;
            }
        }
    }

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_finish);

    for (long long i = i_first; i <= min(i_last, a->last); i++) {
        if (i_primes[i - i_first]) {
            a->prime_count = a->prime_count + 1;
        }
    }
    a->time = double (time_finish.tv_sec - time_start.tv_sec
                     + 1e-9 * (time_finish.tv_nsec - time_start.tv_nsec));
    return a;
}


int main(int argc, char **argv) {
    double max_time = 0, sum_time = 0;
    long long fisrt, last, temp;
    int num_thread, prime_count = 0;
    //ofstream fout(argv[3]);

    pthread_mutex_init(&mutex, NULL);

    sscanf(argv[1], "%llu", &fisrt);
    sscanf(argv[2], "%llu", &last);
    sscanf(argv[3], "%d", &num_thread);

    pthread_t thread_id[num_thread];
    smpl param[num_thread];
    vector<bool> primes((int) sqrt(last + 1.0), true);

    primes[0] = primes[1] = false;

    for (temp = 2; temp * temp <= last; temp++) {
        if (primes[temp]) {
            for (long long j = temp * temp; j * j <= last; j += temp ) {
                primes[j] = false;
            }
        }
    }

    for (long long j = max(fisrt, (long long) 2); j * j <= last; j++) {
        if (primes[j]) {
            prime_count++;
            //fout << j << endl;
        }
    }

    //fout.close();

    for (int i = 0; i < num_thread; i++) {
        param[i] = smpl(num_thread + 1, i + 1, 0, max(fisrt, temp),
                        last, 0, &primes, argv[3]);
        pthread_create(&thread_id[i], NULL, worker, &param[i]);
    }

    for (int i = 0; i < num_thread; i++) {
        smpl *sm_temp;
        pthread_join(thread_id[i], (void **) &sm_temp);
        prime_count += sm_temp->prime_count;
        sum_time += sm_temp->time;
        if (sm_temp->time > max_time) {
        	max_time = sm_temp->time;
        }
    }
    cout << prime_count << endl;
    /*cout << "There are " << prime_count << " primes" << endl
         << "Overall time: " << sum_time << endl
         << "Maximal single process time: " << max_time << endl;*/

    pthread_mutex_destroy(&mutex);

    return 0;
}
