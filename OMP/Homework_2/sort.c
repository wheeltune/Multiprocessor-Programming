#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define swap(type, i, j) {type t = i; i = j; j = t;}

struct InputData {
    int array_size, max_chunk_size, threads_count;
};

int parseInputData(struct InputData* input_data, int argc, char* const argv[]) {
    input_data->array_size = atoi(argv[1]);
    input_data->max_chunk_size = atoi(argv[2]);
    input_data->threads_count = atoi(argv[3]);

    return 0;
}

int qsortCmp (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int upper_bound(const int* begin, const int* end, int target) {
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

void mergeArrays(const int* begin1, const int* end1, const int* begin2, const int* end2, int* destination) {
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

void parallelMerge(const int* begin, const int* middle, const int* end, int* destination) {
    if (end - begin < 2) {
        return;
    }

    int div1 = (middle - begin) >> 1;
    int div2 = upper_bound(middle, end, begin[div1]);

    #pragma omp task
    mergeArrays(begin, begin + div1, middle, middle + div2, destination);

    #pragma omp task
    mergeArrays(begin + div1, middle, middle + div2, end, destination + div1 + div2);
}

void parallelMergeSort(int* begin, int* end, int threads_count, int max_chunk_size) {
    int len = end - begin;
    int step = fmin(max_chunk_size, (int)((len + threads_count - 1) / threads_count));

    int* second_array = (int*) calloc(sizeof(int), len);

    #pragma omp parallel num_threads(threads_count)
    #pragma omp single
    {
        #pragma omp taskgroup
        {
            int i;
            for (i = 0; i < len - step; i += step) {
                #pragma omp task
                qsort(begin + i, step, sizeof(int), qsortCmp);
            }
            if (i < len) {
                #pragma omp task
                qsort(begin + i, len - i, sizeof(int), qsortCmp);
            }
        }

        {
            int i;
            int* source = begin;
            int* destination = second_array;
            for (i = 0; step < len; step <<= 1, ++i) {
                int j;
                #pragma omp taskgroup
                {
                    for (j = 0; j < len - (step << 1); j += (step << 1)) {
                        #pragma omp task
                        parallelMerge(source + j, source + j + step, source + j + (step << 1), destination + j);
                    }
                    if (j + step < len) {
                        #pragma omp task
                        parallelMerge(source + j, source + j + step, source + j + (len - j), destination + j);
                    } else {
                        #pragma omp task
                        memcpy(destination + j, source + j, (len - j) * sizeof(int));
                    }
                }

                swap(int*, source, destination);
            }

            if (i & 1) {
                int j;
                int part_len = (len + threads_count - 1) / threads_count;
                for (j = 0; j < len - part_len; j += part_len) {
                    #pragma omp task
                    memcpy(begin + j, second_array + j, part_len * sizeof(int));
                }
                if (j < len) {
                    #pragma omp task
                    memcpy(begin + j, second_array + j, (len - j) * sizeof(int));
                }
            }
        }
    }

    free(second_array);
}

int main(int argc, char* const argv[]) {
    struct InputData input_data;

    if (parseInputData(&input_data, argc, argv) != 0) {
        return 1;
    }

    int* array  = (int*) calloc(sizeof(int), input_data.array_size);
    int* qarray = (int*) calloc(sizeof(int), input_data.array_size);

    fillRandom(array, input_data.array_size);
    memcpy(qarray, array, input_data.array_size * sizeof(int));

    FILE *data = fopen("data.txt", "a");

    for (size_t i = 0; i < input_data.array_size; ++i) {
        fprintf(data, "%d ", array[i]);
    }
    fprintf(data, "\n");

    // OMP Parralell sort part
    double laps_time = omp_get_wtime();
    parallelMergeSort(array, array + input_data.array_size, input_data.threads_count, input_data.max_chunk_size);
    laps_time = omp_get_wtime() - laps_time;


    FILE* stats = fopen("stats.txt", "a");
    fprintf(stats, "%lfs %d %d %d\n", laps_time,
                                      input_data.array_size,
                                      input_data.max_chunk_size,
                                      input_data.threads_count);
    fclose(stats);

    // QSort part
    double qlaps_time = omp_get_wtime();
    qsort(qarray, input_data.array_size, sizeof(int), qsortCmp);
    qlaps_time = omp_get_wtime() - qlaps_time;


    FILE* qstats = fopen("qstats.txt", "a");
    fprintf(qstats, "%lfs %d\n", qlaps_time,
                                 input_data.array_size);
    fclose(qstats);


    for (size_t i = 0; i < input_data.array_size; ++i) {
        fprintf(data, "%d ", array[i]);
    }
    fprintf(data, "\n\n");

    fclose(data);

    int areEqual = 1;
    for (int i = 0; i < input_data.array_size; ++i) {
        if (array[i] != qarray[i]) {
            areEqual = 0;
            break;
        }
    }
    if (areEqual == 0) {
        printf("WA %d %d %d\n", input_data.array_size, input_data.max_chunk_size, input_data.threads_count);
    }

    free(array);
    free(qarray);

    return 0;
}
