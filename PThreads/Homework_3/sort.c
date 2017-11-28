#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "threadpool.h"
#include "log.h"

#define swap(type, i, j) {type t = i; i = j; j = t;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct InputData {
    int array_size, max_chunk_size, threads_count;
};

int parseInputData(struct InputData* input_data, int argc, char* const argv[]) {
    input_data->array_size = atoi(argv[1]);
    input_data->max_chunk_size = atoi(argv[2]);
    input_data->threads_count = atoi(argv[3]);

    return 0;
}

void get_time(struct timespec *to_store) {
    clock_gettime(CLOCK_MONOTONIC, to_store);
}

double elapsed_time(struct timespec from) {
    struct timespec now;
    get_time(&now);

    double elapsed = (now.tv_sec - from.tv_sec);
    elapsed += (now.tv_nsec - from.tv_nsec) / 1000000000.0;
    return elapsed;
}

int qsortCmp (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int upperBound(const int* begin, const int* end, int target) {
    int left = 0, right = end - begin;
    while (right - left > 1) {
        int middle = (right + left) >> 1;
        if (begin[middle] > target) {
            right = middle;
        } else {
            left = middle;
        }
    }
    if (begin[left] > target) {
        return left;
    } else {
        return right;
    }
}

void fillRandom(int* array, int len) {
    srand(time(NULL));
    for (int i = 0; i < len; ++i) {
        array[i] = rand() % 100;
    }
}

// MERGE ARRAYS
void mergeArrays(const int* begin1,
                 const int* end1,
                 const int* begin2,
                 const int* end2,
                 int* destination) {
    const int* it1 = begin1;
    const int* it2 = begin2;

    while(it1 != end1 && it2 != end2) {
        if (*it1 < *it2) {
            *destination = *it1;
            it1++;
        } else {
            *destination = *it2;
            it2++;
        }
        destination++;
    }

    while (it1 != end1) {
        *destination = *it1;
        it1++;
        destination++;
    }

    while (it2 != end2) {
        *destination = *it2;
        it2++;
        destination++;
    }
}

struct ThreadMergeArraysArg {
    const int *begin1;
    const int *end1;
    const int *begin2;
    const int *end2;
    int *destination;
};

void threadMergeArrays(void *arg) {
    struct ThreadMergeArraysArg *args = (struct ThreadMergeArraysArg *)arg;
    mergeArrays(args->begin1,
                args->end1,
                args->begin2,
                args->end2,
                args->destination);
    free(args);
}

// PARALLEL MERGE
void parallelMerge(struct ThreadPool *thread_pool,
                   const int* begin,
                   const int* middle,
                   const int* end,
                   int* destination) {
    if (end - begin < 2) {
        return;
    }

    int div1 = (middle - begin) >> 1;
    int div2 = upperBound(middle, end, begin[div1]);

    struct ThreadMergeArraysArg *args =
        (struct ThreadMergeArraysArg*) malloc(sizeof(struct ThreadMergeArraysArg));
    args->begin1      = begin;
    args->end1        = begin + div1;
    args->begin2      = middle;
    args->end2        = middle + div2;
    args->destination = destination;
    ThreadPool_add(thread_pool, threadMergeArrays, (void*)args);

    args =
        (struct ThreadMergeArraysArg*) malloc(sizeof(struct ThreadMergeArraysArg));
    args->begin1      = begin + div1;
    args->end1        = middle;
    args->begin2      = middle + div2;
    args->end2        = end;
    args->destination = destination + div1 + div2;
    ThreadPool_add(thread_pool, threadMergeArrays, (void*)args);
}

struct ThreadParallelMergeArg {
    struct ThreadPool *thread_pool;
    int *begin;
    int *middle;
    int *end;
    int *destination;
};

void threadParallelMerge(void *arg) {
    struct ThreadParallelMergeArg *args = (struct ThreadParallelMergeArg*) arg;
    parallelMerge(args->thread_pool,
                  args->begin,
                  args->middle,
                  args->end,
                  args->destination);
    free(args);
}

// QSORT
struct ThreadQSortArg {
    void* base;
    size_t num;
    size_t size;
    int (*compar)(const void*, const void*);
};

void threadQSort(void *arg) {
    struct ThreadQSortArg *args = (struct ThreadQSortArg *)arg;
    qsort(args->base,
          args->num,
          args->size,
          args->compar);
    free(args);
}

// MEMCPY
struct ThreadMemCpyArg {
    void *str1;
    const void *str2;
    size_t n;
};

void threadMemCpy(void *arg) {
    struct ThreadMemCpyArg *args = (struct ThreadMemCpyArg*) arg;
    memcpy(args->str1, args->str2, args->n);
    free(args);
}


void parallelMergeSort(int* begin, int* end, int threads_count, int max_chunk_size) {
    struct ThreadPool thread_pool;
    ThreadPool_construct(&thread_pool, threads_count);

    int len = end - begin;
    int step = MIN(max_chunk_size, (int)((len + threads_count - 1) / threads_count));

    int* second_array = (int*) calloc(sizeof(int), len);

    {
        int i;
        for (i = 0; i < len - step; i += step) {
            struct ThreadQSortArg *args =
                (struct ThreadQSortArg*) malloc(sizeof(struct ThreadQSortArg));
            args->base   = begin + i;
            args->num    = step;
            args->size   = sizeof(int);
            args->compar = qsortCmp;

            ThreadPool_add(&thread_pool, threadQSort, (void*)args);
        }
        if (i < len) {
            struct ThreadQSortArg *args =
                (struct ThreadQSortArg*) malloc(sizeof(struct ThreadQSortArg));
            args->base   = begin + i;
            args->num    = len - i;
            args->size   = sizeof(int);
            args->compar = qsortCmp;

            ThreadPool_add(&thread_pool, threadQSort, (void*)args);
        }
    }

    printf("Before wait");
    fflush(stdout);
    ThreadPool_wait(&thread_pool);
    printf("After wait");
    fflush(stdout);

    {
        int i;
        int* source = begin;
        int* destination = second_array;
        for (i = 0; step < len; step <<= 1, ++i) {
            int j;
            for (j = 0; j < len - (step << 1); j += (step << 1)) {
                struct ThreadParallelMergeArg *args =
                    (struct ThreadParallelMergeArg*) malloc(sizeof(struct ThreadParallelMergeArg));
                args->thread_pool = &thread_pool;
                args->begin       = source + j;
                args->middle      = source + j + step;
                args->end         = source + j + (step << 1);
                args->destination = destination + j;

                ThreadPool_add(&thread_pool, threadParallelMerge, (void*)args);
            }
            if (j + step < len) {
                struct ThreadParallelMergeArg *args =
                    (struct ThreadParallelMergeArg*) malloc(sizeof(struct ThreadParallelMergeArg));
                args->thread_pool = &thread_pool;
                args->begin       = source + j;
                args->middle      = source + j + step;
                args->end         = source + j + (len - j);
                args->destination = destination + j;

                ThreadPool_add(&thread_pool, threadParallelMerge, (void*)args);
            } else {
                struct ThreadMemCpyArg *args =
                    (struct ThreadMemCpyArg*) malloc(sizeof(struct ThreadMemCpyArg));
                args->str1 = destination + j;
                args->str2 = source + j;
                args->n    = (len - j) * sizeof(int);

                ThreadPool_add(&thread_pool, threadMemCpy, (void*) args);
            }
            ThreadPool_wait(&thread_pool);

            swap(int*, source, destination);
        }

        if (i & 1) {
            int j;
            int part_len = (len + threads_count - 1) / threads_count;
            for (j = 0; j < len - part_len; j += part_len) {
                struct ThreadMemCpyArg *args =
                    (struct ThreadMemCpyArg*) malloc(sizeof(struct ThreadMemCpyArg));
                args->str1 = begin + j;
                args->str2 = second_array + j;
                args->n    = part_len * sizeof(int);

                ThreadPool_add(&thread_pool, threadMemCpy, (void*) args);
            }
            if (j < len) {
                struct ThreadMemCpyArg *args =
                    (struct ThreadMemCpyArg*) malloc(sizeof(struct ThreadMemCpyArg));
                args->str1 = begin + j;
                args->str2 = second_array + j;
                args->n    = (len - j) * sizeof(int);

                ThreadPool_add(&thread_pool, threadMemCpy, (void*) args);
            }
        }
    }
    ThreadPool_join(&thread_pool);
    free(second_array);

    ThreadPool_destruct(&thread_pool);
}

int main(int argc, char* const argv[]) {
    struct InputData input_data;

    ENABLE_LOGGING = 1;

    if (parseInputData(&input_data, argc, argv) != 0) {
        return 1;
    }

    int* array  = (int*) calloc(sizeof(int), input_data.array_size);
    int* qarray = (int*) calloc(sizeof(int), input_data.array_size);
    int* tmp    = (int*) calloc(sizeof(int), input_data.array_size);

    fillRandom(array, input_data.array_size);
    memcpy(qarray, array, input_data.array_size * sizeof(int));
    memcpy(tmp, array, input_data.array_size * sizeof(int));
    // Parralell sort part
    printf("Start");
    fflush(stdout);
    struct timespec start_time;
    get_time(&start_time);
    parallelMergeSort(array,
                      array + input_data.array_size,
                      input_data.threads_count,
                      input_data.max_chunk_size);
    double laps_time = elapsed_time(start_time);
    printf("End");
    fflush(stdout);

    FILE* stats = fopen("stats.txt", "a");
    fprintf(stats, "%lfs %d %d %d\n", laps_time,
                                      input_data.array_size,
                                      input_data.max_chunk_size,
                                      input_data.threads_count);
    fclose(stats);

    // QSort part
    struct timespec qstart_time;
    get_time(&qstart_time);
    qsort(qarray, input_data.array_size, sizeof(int), qsortCmp);
    double qlaps_time = elapsed_time(qstart_time);

    FILE* qstats = fopen("qstats.txt", "a");
    fprintf(qstats, "%lfs %d\n", qlaps_time,
                                 input_data.array_size);
    fclose(qstats);

    int areEqual = 1;
    for (int i = 0; i < input_data.array_size; ++i) {
        if (array[i] != qarray[i]) {
            areEqual = 0;
            break;
        }
    }
    if (areEqual == 0) {
        printf("WA %d %d %d\n", input_data.array_size,
                                input_data.max_chunk_size,
                                input_data.threads_count);
    }

    free(array);
    free(qarray);
    free(tmp);

    return 0;
}
