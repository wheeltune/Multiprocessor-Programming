#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum GeneratorState { STOP, RUN } generator_state_ = STOP;

const size_t CHUNKS_COUNT = 5000;
const size_t CHUNK_SIZE = 5000;

double** chunks_;
size_t chunk_num_ = 0;

long long summary_time_ = 0;
int finised_in_b_ = 0;
double laps_time_ = 0.0;

int a_, b_, x_, N_;
size_t num_threads_;
double p_;

double randomDouble() {
    return (double)rand() / (double)RAND_MAX;
}

double* generateNewChunk() {
    double* new_chunk = (double*) malloc(sizeof(double) * CHUNK_SIZE);
    for (size_t i = 0; i < CHUNK_SIZE; ++i) {
        new_chunk[i] = randomDouble();
    }
    return new_chunk;
}

void constructGenerator() {
    chunks_ = (double**) malloc(sizeof(double*) * CHUNKS_COUNT);
    for (size_t i = 0; i < CHUNKS_COUNT; ++i) {
        chunks_[i] = generateNewChunk();
    }
}

void destructGenerator() {
    for (size_t i = 0; i < CHUNKS_COUNT; ++i) {
        if (chunks_[i] != NULL) {
            free(chunks_[i]);
        }
    }
    free(chunks_);
}

void runGenerator() {
    generator_state_ = RUN;

    srand(time(NULL));
    while (generator_state_ == RUN) {
        for (size_t i = 0; i < CHUNKS_COUNT && generator_state_ == RUN; ++i) {
            // if (chunks_[i] == NULL) {
            //     chunks_[i] = generateNewChunk();
            // }

            while (chunks_[i] != NULL && generator_state_ == RUN) {}

            if (generator_state_ == STOP) {
                break;
            }
            chunks_[i] = generateNewChunk();
        }
    }
}

void stopGenerator() {
    generator_state_ = STOP;
}

double* getNewChunk() {
    if (num_threads_ == 1) {
        return generateNewChunk();
    }
    
    double* new_chunk = NULL;
    #pragma omp critical
    {
        while (chunks_[chunk_num_] == NULL) {}

        new_chunk = chunks_[chunk_num_];
        chunks_[chunk_num_] = NULL;

        chunk_num_++;
        if (chunk_num_ >= CHUNKS_COUNT) {
            chunk_num_ -= CHUNKS_COUNT;
        }
    }
    return new_chunk;
}

double getRandom(size_t* pos, double** local_chunk) {
    if (*pos >= CHUNK_SIZE) {
        free(*local_chunk);
        *local_chunk = getNewChunk();
        *pos = 0;
    }
    return (*local_chunk)[(*pos)++];
}

int checkParameters() {
    if (a_ >= b_) {
        printf("'a' must be grater than 'b'\n");
        return 1;
    }
    if (x_ <= a_ || x_ >= b_) {
        printf("'x' should be inside ('a', 'b')\n");
        return 1;
    }
    if (N_ <= 0) {
        printf("'N' should be positive\n");
        return 1;
    }
    if (p_ < 0.0 || p_ > 1.0) {
        printf("'x' should be inside [0, 1]\n");
        return 1;
    }
    return 0;
}

int main(int argc, char* argv[]) {
    a_ = atoi(argv[1]); b_ = atoi(argv[2]);
    x_ = atoi(argv[3]); N_ = atoi(argv[4]);
    p_ = atof(argv[5]); num_threads_ = atoi(argv[6]);

    if (checkParameters() > 0) {
        return 1;
    }

    omp_set_nested(1);

    int thread_team1 = 0;
    int thread_team2 = num_threads_;
    if (thread_team2 > 1) {
        thread_team1++;
        thread_team2--;
    }

    constructGenerator();

    #pragma omp parallel num_threads(thread_team1 + 1) 
    {
        if (omp_get_thread_num() == 1) {
            #pragma omp parallel num_threads(thread_team1)
            {
                runGenerator();
            }
        } else {
            double start_time = omp_get_wtime();
            #pragma omp parallel num_threads(thread_team2) reduction(+:finised_in_b_) reduction(+:summary_time_)
            {
                double* local_chunk = getNewChunk();
                size_t chunk_pos = 0;

                #pragma omp for schedule(dynamic)
                for (int i = 0; i < N_; ++i) {                  
                    int current_pos = x_;

                    while (current_pos != a_ && current_pos != b_) {
                        if (getRandom(&chunk_pos, &local_chunk) < p_) {
                            current_pos++;
                        } else {
                            current_pos--;
                        }
                        summary_time_++;
                    }

                    if (current_pos == b_) {
                        finised_in_b_++;
                    }
                }

                free(local_chunk);   
            }
            laps_time_ = omp_get_wtime() - start_time;
            stopGenerator();
        }
    }

    destructGenerator();

    FILE* stats = fopen("stats.txt", "a");
    fprintf(stats, "%lf %lf %lfs %d %d %d %d %lf %zu\n", (double)finised_in_b_ / (double)N_ * 100, 
                                                         (double)summary_time_ / (double)N_,
                                                         laps_time_,
                                                         a_, b_, x_, N_, p_, num_threads_);
    fclose(stats);

    printf("%lf %lf %lfs %d %d %d %d %lf %zu\n", (double)finised_in_b_ / (double)N_ * 100, 
                                                 (double)summary_time_ / (double)N_,
                                                 laps_time_,
                                                 a_, b_, x_, N_, p_, num_threads_);
    return 0;
}