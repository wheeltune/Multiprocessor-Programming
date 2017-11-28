#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MASTER 0
#define CHUNK_LENGTH 500
#define CHUNKS_COUNT 500

enum side_e {
    LEFT = 0, RIGHT = 1, UP = 2, DOWN = 3, NONE = 4
};

struct point_t {
    int x, y;
};

struct chunk_t {
    int* data;
    int pos;
};

void chunk_init(struct chunk_t* chunk) {
    chunk->data = malloc(CHUNK_LENGTH * sizeof(int));
    chunk->pos = CHUNK_LENGTH;
}

void chunk_destroy(struct chunk_t* chunk) {
    free(chunk->data);
}

MPI_Datatype point_mpi_type;
void point_mpi_type_init(MPI_Datatype* type) {
    int blocks[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Aint displacements[2];

    MPI_Aint intex;
    MPI_Type_extent(MPI_INT, &intex);
    displacements[0] = (MPI_Aint) 0;
    displacements[1] = intex;

    MPI_Type_struct(2, blocks, displacements, types, type);
    MPI_Type_commit(type);
}

int l_, a_, b_, n_, N_;
double p_l_, p_r_, p_u_, p_d_;

double start_time;

void run_generator(int count_fields) {
    int** chunks_data = (int**) calloc(CHUNKS_COUNT, sizeof(int*));
    int count_finished = 0;
    #pragma omp parallel shared(count_finished, chunks_data) num_threads(2)
    {
        if (omp_get_thread_num() == 0) {
            srand(time(NULL));
            int* finished = (int*) calloc(count_fields, sizeof(int));
            int chunk_pos = 0;

            int interlocutor_rank = 1;
            while (count_finished < count_fields) {
                if (!finished[interlocutor_rank - 1]) {
                    int flag = 0;
                    MPI_Status status;
                    MPI_Iprobe(interlocutor_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
                    if (flag) {
                        int is_finished = 0;
                        MPI_Recv(&is_finished, 1, MPI_INT, interlocutor_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        if (is_finished) {
                            #pragma omp atomic
                            count_finished += 1;

                            finished[interlocutor_rank - 1] = 1;
                            continue;
                        }

                        while (chunks_data[chunk_pos] == NULL) {}
                        MPI_Ssend(chunks_data[chunk_pos], CHUNK_LENGTH, MPI_INT, interlocutor_rank, 0, MPI_COMM_WORLD);
                        free(chunks_data[chunk_pos]);
                        chunks_data[chunk_pos] = NULL;

                        chunk_pos++;
                        if (chunk_pos == CHUNKS_COUNT) {
                            chunk_pos = 0;
                        }
                    }
                }
                interlocutor_rank++;
                if (interlocutor_rank > count_fields) {
                    interlocutor_rank = 1;
                }
            }

            free(finished);

         } else {
            int chunk_pos = 0;
            while (count_finished < count_fields) {
                int* new_chunk = (int*) malloc(CHUNK_LENGTH * sizeof(int));
                for (int i = 0; i < CHUNK_LENGTH; ++i) {
                    new_chunk[i] = rand();
                }

                while (chunks_data[chunk_pos] != NULL && count_finished < count_fields) {}
                chunks_data[chunk_pos] = new_chunk;

                chunk_pos++;
                if (chunk_pos >= CHUNKS_COUNT) {
                    chunk_pos = 0;
                }
            }
        }
    }

    for (int i = 0; i < CHUNKS_COUNT; ++i) {
        if (chunks_data[i] != NULL) {
            free(chunks_data[i]);
        }
    }
    free(chunks_data);
    double laps_time = MPI_Wtime() - start_time;

    FILE* stats_file = fopen("stats.txt", "w");
    fprintf(stats_file, "%d %d %d %d %d %lf %lf %lf %lf %lfs\n", l_, a_, b_, n_, N_, p_l_, p_r_, p_u_, p_d_, laps_time);
    for (int i = 0; i < count_fields; ++i) {
        int count_points;
        MPI_Recv(&count_points, 1, MPI_INT, i + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        fprintf(stats_file, "%d: %d\n", i, count_points);
    }
    fclose(stats_file);
}

int get_random(struct chunk_t* chunk) {
    if (chunk->pos >= CHUNK_LENGTH) {
        int is_finished = 0;
        MPI_Ssend(&is_finished, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

        MPI_Status status;
        MPI_Recv(chunk->data, CHUNK_LENGTH, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        chunk->pos = 0;
    }
    return chunk->data[chunk->pos++];
}

double get_random_double(struct chunk_t* chunk) {
    return (double) get_random(chunk) / (double) RAND_MAX;
}

void run_process(int rank, int count_fields) {
    int ranks[4];
    ranks[LEFT] = rank - 1;
    if (ranks[LEFT] % a_ == 0) {
        ranks[LEFT] += a_;
    }
    ranks[RIGHT] = rank + 1;
    if (ranks[RIGHT] % a_ == 1) {
        ranks[RIGHT] -= a_;
    }
    ranks[UP] = rank + a_;
    if (ranks[UP] > count_fields) {
        ranks[UP] -= count_fields;
    }
    ranks[DOWN] = rank - a_;
    if (ranks[DOWN] <= 0) {
        ranks[DOWN] += count_fields;
    }

    struct point_t min_point, max_point;
    {
        int horizontal_index = (rank - 1) % a_;
        int vertical_index = (rank - 1) / a_;
        min_point.x = horizontal_index * l_;
        min_point.y = vertical_index * l_;
        max_point.x = min_point.x + l_;
        max_point.y = min_point.y + l_;
    }

    struct chunk_t chunk;
    chunk_init(&chunk);
    struct point_t* points = (struct point_t*) malloc(N_ * sizeof(struct point_t));
    int count_points = N_;
    for (int i = 0; i < N_; ++i) {
        points[i].x = get_random(&chunk) % l_ + min_point.x;
        points[i].y = get_random(&chunk) % l_ + min_point.y;
    }

    for (int iter = 0; iter < n_; ++iter) {
        struct point_t* to[5];
        int to_count[5] = {0, 0, 0, 0, 0};
        for (int side = 0; side < 4; ++side) {
            to[side]  = (struct point_t*) malloc(count_points * sizeof(struct point_t));
        }
        to[NONE]  = (struct point_t*) malloc(N_ * count_fields * sizeof(struct point_t));

        for (int i = 0; i < count_points; ++i) {
            double random = get_random_double(&chunk);
            if (random <= p_l_) {
                points[i].x--;
                if (points[i].x < 0) {
                    points[i].x = l_ * a_ - 1;
                }
                if (points[i].x < min_point.x || points[i].x >= max_point.x) {
                    memcpy(to[LEFT] + to_count[LEFT], points + i, sizeof(struct point_t));
                    to_count[LEFT]++;
                    continue;
                }
            } else if (random <= p_l_ + p_r_) {
                points[i].x++;
                if (points[i].x >= l_ * a_) {
                    points[i].x = 0;
                }
                if (points[i].x < min_point.x || points[i].x >= max_point.x) {
                    memcpy(to[RIGHT] + to_count[RIGHT], points + i, sizeof(struct point_t));
                    to_count[RIGHT]++;
                    continue;
                }
            } else if (random <= p_l_ + p_r_ + p_d_) {
                points[i].y--;
                if (points[i].y < 0) {
                    points[i].y = l_ * b_ - 1;
                }
                if (points[i].y < min_point.y || points[i].y >= max_point.y) {
                    memcpy(to[DOWN] + to_count[DOWN], points + i, sizeof(struct point_t));
                    to_count[DOWN]++;
                    continue;
                }
            } else {
                points[i].y++;
                if (points[i].y >= l_ * b_) {
                    points[i].y = 0;
                }
                if (points[i].y < min_point.y || points[i].y >= max_point.y) {
                    memcpy(to[UP] + to_count[UP], points + i, sizeof(struct point_t));
                    to_count[UP]++;
                    continue;
                }
            }
            memcpy(to[NONE] + to_count[NONE], points + i, sizeof(struct point_t));
            to_count[NONE]++;
        }

        MPI_Request requests[4];
        for (int side = 0; side < 4; ++side) {
            MPI_Issend(to[side], to_count[side], point_mpi_type, ranks[side], 0, MPI_COMM_WORLD, requests + side);
        }

        MPI_Status status;
        for (int side = 0; side < 4; ++side) {
            MPI_Probe(ranks[side], MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            int recv_count;
            MPI_Get_count(&status, point_mpi_type, &recv_count);

            struct point_t* recv_points = (struct point_t*) malloc(recv_count * sizeof(struct point_t));
            MPI_Recv(recv_points, recv_count, point_mpi_type, ranks[side], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            memcpy(to[NONE] + to_count[NONE], recv_points, recv_count * sizeof(struct point_t));
            to_count[NONE] += recv_count;
        }
        MPI_Waitall(4, requests, MPI_STATUS_IGNORE);

        free(points);
        points = to[NONE];
        count_points = to_count[NONE];
    }

    {
        int is_finished = 1;
        MPI_Ssend(&is_finished, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
        MPI_Ssend(&count_points, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
    }
    free(points);
    chunk_destroy(&chunk);
}

int main(int argc, char **argv) {
    l_ = atoi(argv[1]);
    a_ = atoi(argv[2]), b_ = atoi(argv[3]);
    n_ = atoi(argv[4]), N_ = atoi(argv[5]);
    p_l_ = atof(argv[6]), p_r_ = atof(argv[7]);
    p_u_ = atof(argv[8]), p_d_ = atof(argv[9]);

    MPI_Init(&argc, &argv);
    point_mpi_type_init(&point_mpi_type);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start_time = MPI_Wtime();
    if (rank == MASTER) {
        run_generator(size - 1);
    } else {
        run_process(rank, size - 1);
    }

    MPI_Finalize();
    return 0;
}
