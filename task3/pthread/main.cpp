
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <time.h>


using namespace std;

pthread_mutex_t mutex;

struct smpl{
    long long start;
    long long end;
    int size;
    int rank;
    int count;
    double time;
    std::vector<bool> * prime;
    char *file;
    smpl(int s = 1 , int r = 1, int c = 0, long long st = 0,
         long long e = 1, double t = 0, std::vector<bool> *p = NULL,
         char *f = 0) {
        size = s;
        rank = r;
        count = c;
        start = st;
        end = e;
        time = t;
        prime = p;
        file = f;
    }
};

void *worker(void *atr) {

    smpl *a = (smpl*) atr;
    long long ibegin = (long long) (a->rank - 1) * (a->end - a->start + 1)
                        / (a->size - 1) + a->start;
    long long iend = (long long ) a->rank * (a->end - a->start + 1)
                      / (a->size - 1) + a->start- 1;
    std::vector<bool> prime1(iend - ibegin + 1, true);
    timespec time_start;
    timespec time_finish;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_start);

    for (long long j = 2; j * j <= a->end; j++) {
        if((*(*a).prime)[j]) {
            for (long long i = (ibegin / j + 1 * (ibegin % j != 0)) * j;
                i <= min(iend, a->end); i += j) {
                prime1[i - ibegin] = false;
            }
        }

    }

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_finish);
    for (long long i = ibegin; i <= min(iend, a->end); i++) {

        if (prime1[i - ibegin]) {
            a->count = a->count + 1;
        }
    }
    a->time = double (time_finish.tv_sec - time_start.tv_sec
                     + 1e-9 * (time_finish.tv_nsec - time_start.tv_nsec));
    return a;
}


int main(int argc, char **argv) {

    long long begin, end;
    pthread_mutex_init(&mutex, NULL);
    sscanf(argv[1], "%llu", &begin);
    sscanf(argv[2], "%llu", &end);
    int numthread;
    sscanf(argv[4], "%d", &numthread);
    std::vector<bool> prime((int) sqrt(end + 1.0), true);
    prime[0] = prime[1] = false;
    long long temp;
    int count = 0;
    ofstream fout(argv[3]);
    for (temp = 2; temp * temp <= end; temp++) {
        if (prime[temp]) {
            for (long long j = temp * temp; j * j <= end; j += temp ) {
                prime[j] = false;
            }
        }

    }
    for (long long j = max(begin, (long long) 2); j * j <= end; j++) {
        if (prime[j]) {
            count++;
            fout << j << endl;
        }
    }
    fout.close();
    pthread_t tid[numthread];
    smpl param[numthread];
    for (int i = 0; i < numthread; i++) {
        param[i] = smpl(numthread + 1, i + 1, 0, max(begin, temp), end, 0, &prime, argv[3]);
        pthread_create(&tid[i], NULL, worker, &param[i]);
    }
    double maxtime = 0;
    double sumtime = 0;
    for (int i = 0; i < numthread; i++) {
        smpl *smtemp;
        pthread_join(tid[i],(void **) &smtemp);
        count += smtemp->count;
        sumtime += smtemp->time;
        if (smtemp->time > maxtime) {
        	maxtime = smtemp->time;
        }
    }
    cout << count << endl;
    cout << "maxtime: " << maxtime << endl << "sumtime: " << sumtime << endl;
    pthread_mutex_destroy(&mutex);
    return 0;
}
